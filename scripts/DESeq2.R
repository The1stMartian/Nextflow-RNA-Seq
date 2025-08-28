# Run DESeq2 on a counts matrix
# Inputs are (2):
# 1) metadata file with grouping information
# 2) a single counts matrix file

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(DESeq2)
  library(ggplot2)
  library(stringr)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

###########################################
# Significance cutoff
padj_cutoff = 0.05

# Volcano Plot Threshold Cutoffs
lfc_thr  <- 1         # |log2FC| cutoff
padj_thr <- 0.05      # FDR cutoff
volcanoPlotName = "VolcanoPlot.png"
###########################################

# --- Ingest counts and metadata
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript runDESeq2.R <counts_matrix.csv> <metadata.csv>")
}

counts_path <- args[1]
meta_path   <- args[2]

# --- Load featureCounts matrix ---
# Handles standard featureCounts columns: Geneid Chr Start End Strand Length + sample columns
message("Reading counts: ", counts_path)
counts_df <- readr::read_csv(counts_path, comment = "#", show_col_types = FALSE)

# Identify sample count columns (everything except the 6 standard columns if present)
std_cols <- intersect(colnames(counts_df), c("Geneid","Chr","Start","End","Strand","Length"))
sample_cols <- setdiff(colnames(counts_df), std_cols)
if (length(sample_cols) < 2) stop("No sample columns detected in the featureCounts file.")

# Tidy column names to sample basenames (strip paths and .bam)
tidy_name <- function(x) {
  x <- sub("^.*/", "", x)       # drop directory
  x <- sub("\\.bam$", "", x)    # drop .bam
  x
}

clean_sample_names <- vapply(sample_cols, tidy_name, character(1))
colnames(counts_df)[match(sample_cols, colnames(counts_df))] <- clean_sample_names

# Build count matrix
gene_col <- if ("Geneid" %in% colnames(counts_df)) "Geneid" else colnames(counts_df)[1]
mat <- as.matrix(counts_df[, clean_sample_names, drop = FALSE])
rownames(mat) <- counts_df[[gene_col]]
storage.mode(mat) <- "integer"

#######################################################
# --- Load metadata ---
message("Reading metadata: ", meta_path)
meta_guess <- tools::file_ext(meta_path)
if (tolower(meta_guess) %in% c("csv")) {
  meta <- readr::read_csv(meta_path, show_col_types = FALSE)
} else {
  meta <- readr::read_tsv(meta_path, show_col_types = FALSE)
}

derive_sample <- function(fq_path) {
  b <- basename(fq_path)
  # strip compression + fastq/fq
  b <- sub("\\.(fastq|fq)(\\.(gz|bz2|xz))?$", "", b, ignore.case = TRUE)
  # remove common Illumina suffix blocks from the end:
  # e.g., _S1_L001_R1_001  or _L001_R1_001  or _R1_001  or _R1  or _1
  b <- sub("(_S\\d+)?(_L\\d{3})?(_R?[12]|_[12])(_\\d{3})?$", "", b, perl = TRUE)
  b
}

meta <- tibble(
  sample    = vapply(meta$fastq_1, derive_sample, character(1)),
  condition = tolower(as.character(meta$condition))
)

required_cols <- c("sample", "condition")
if (!all(required_cols %in% colnames(meta))) {
  stop("Metadata must contain columns: 'sample' and 'condition'")
}

meta <- meta %>%
  mutate(sample = as.character(sample),
         condition = as.character(condition))

# Keep only samples present in counts and in metadata (and order them)
present <- meta$sample %in% colnames(mat)
if (!all(present)) {
  missing <- meta$sample[!present]
  stop(paste0("Samples in metadata not found in counts columns: ", paste(missing, collapse = ", ")))
}
mat <- mat[, meta$sample, drop = FALSE]

# Condition factor with explicit control->test order
if (!all(meta$condition %in% c("control","test"))) {
  stop("Metadata 'condition' must be either 'control' or 'test'.")
}
meta$condition <- factor(meta$condition, levels = c("control","test"))
rownames(meta) <- meta$sample

# --- Gene filtering: remove very low counts ---
keep_genes <- rowSums(mat) >= 10
mat_f <- mat[keep_genes, , drop = FALSE]
if (nrow(mat_f) == 0) stop("No genes pass the minimal count filter (rowSums >= 10).")


# --- Gene filtering: remove very low counts ---
keep_genes <- rowSums(mat) >= 10
mat_f <- mat[keep_genes, , drop = FALSE]
if (nrow(mat_f) == 0) stop("No genes pass the minimal count filter (rowSums >= 10).")

# --- DESeq2 ---
dds <- DESeqDataSetFromMatrix(countData = mat_f,
                              colData   = meta,
                              design    = ~ condition)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "test", "control"))
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  arrange(padj, pvalue)

# Write results of DE
readr::write_tsv(res_df, "deseq2_results_test_vs_control.tsv")

# Significant genes (padj < 0.05 - adjustable at the top)
signif_df <- res_df %>%
  filter(!is.na(padj), padj < padj_cutoff)
readr::write_tsv(signif_df, "deseq2_results_test_vs_control.signif.tsv")

# Normalized counts and size factors
norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>% rownames_to_column("gene")
readr::write_tsv(norm_counts, "deseq2_normalized_counts.tsv")

sf <- sizeFactors(dds)
readr::write_tsv(tibble(sample = names(sf), size_factor = as.numeric(sf)),
                 "size_factors.tsv")

# Plots: dispersion, MA, PCA (rlog)
rld <- rlog(dds, blind = TRUE)
pdf("deseq2_plots.pdf")
plotDispEsts(dds)
plotMA(dds, main = "DESeq2 MA: test vs control")
print(plotPCA(rld, intgroup = "condition"))
dev.off()

# Session info
writeLines(c(capture.output(sessionInfo())), "deseq2_sessionInfo.txt")

########################################################################
# Volcano plot

# ---- identify columns ----
res = res_df
gene_col <- c("gene", "Geneid", "Gene", "symbol")
gene_col <- gene_col[gene_col %in% names(res)][1]

if (is.na(gene_col)) stop("Could not find a gene column (tried: gene, Geneid, Gene, symbol).")

if (!"log2FoldChange" %in% names(res)) stop("Missing 'log2FoldChange' column.")
if (!"padj" %in% names(res) && !"pvalue" %in% names(res)) {
  stop("Need at least 'padj' or 'pvalue' column.")
}
if (!"padj" %in% names(res)) {
  warning("No 'padj' column found; using 'pvalue' instead (no multiple-testing correction).")
  res$padj <- res$pvalue
}

# ---- prep data ----
df <- res %>%
  transmute(
    gene = .data[[gene_col]],
    log2FC = as.numeric(log2FoldChange),
    padj   = as.numeric(padj),
    neglog10padj = -log10(padj)
  )

# Log2FC and FDR set above
df <- df %>%
  mutate(
    sig = case_when(
      !is.na(padj) & padj <= padj_thr & log2FC >=  lfc_thr ~ "Up",
      !is.na(padj) & padj <= padj_thr & log2FC <= -lfc_thr ~ "Down",
      TRUE ~ "NS"
    )
  )

# axis limits (symmetric around 0, robust to outliers)
highest_log2FC <- max(abs(df$log2FC[is.finite(df$log2FC)]))
xlim_val <- highest_log2FC +0.5

# choose genes to label: top by smallest padj among significant
top_n_label  <- 15 # label the top X most significant genes

to_label <- df %>%
  filter(sig %in% c("Up","Down"), is.finite(neglog10padj)) %>%
  arrange(padj, desc(abs(log2FC))) %>%
  slice_head(n = top_n_label)

maxY = 1+ max(df$neglog10padj)

# ---- plot ----
p <- ggplot(df, aes(x = log2FC, y = neglog10padj)) +
  geom_point(aes(shape = sig, alpha = sig), size = 2.5) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2, linewidth = 0.5) +
  geom_hline(yintercept = -log10(padj_thr), linetype = 2, linewidth = 0.5) +
  geom_text_repel(
    data = to_label,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    min.segment.length = 0
  ) +
  scale_shape_manual(values = c(NS = 16, Up = 17, Down = 25)) +
  scale_alpha_manual(values = c(NS = 0.3, Up = 0.9, Down = 0.9), guide = "none") +
  coord_cartesian(xlim = c(-xlim_val, xlim_val), ylim = c(0, maxY)) +
  labs(
    title = "Volcano plot",
    subtitle = sprintf("padj ≤ %.3g, |log2FC| ≥ %.3g", padj_thr, lfc_thr),
    x = "log2(Fold Change)",
    y = expression(-log[10]("padj")),
    shape = NULL
  ) +
  theme_minimal(base_size = 12) +
  scale_shape_manual(values = c(NS = 16, Up = 17, Down = 25)) +
  scale_alpha_manual(values = c(NS = 0.3, Up = 0.9, Down = 0.9), guide = "none") +
  scale_color_manual(values = c(NS = "grey70", Up = "blue", Down = "red")) 

# ---- save PNG & PDF ----
ggsave(volcanoPlotName, plot = p, width = 7, height = 5, dpi = 300)
ggsave(sub("\\.png$", ".pdf", volcanoPlotName, ignore.case = TRUE), plot = p, width = 7, height = 5)

message("Saved: ", volcanoPlotName)