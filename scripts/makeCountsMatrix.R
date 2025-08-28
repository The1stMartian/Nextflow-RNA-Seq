# Summary: From counts file list as args, outputs a combined counts matrix
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
})

# Default
outputCountsFile <- "countsMatrix.csv"

# Grab file list and output file name from command line
args <- commandArgs(trailingOnly = TRUE)
outputCountsFile <- args[length(args)]        # output file name
args <- args[-length(args)]                   # drop last entry (outfile name)

if (length(args) < 2) {
  stop("Usage: Rscript combine_featurecounts.R file1.featurecounts.txt file2.featurecounts.txt ...")
}

# Gene file names from paths 
fn = sapply(args, basename, USE.NAMES = FALSE)

# Helper: extract sample name from filename (e.g. test3.featurecounts.txt -> test3)
get_sample_name <- function(fname) {
  name =  unlist(strsplit(fname, ".", fixed=T))[[1]]
  return(name)
}

sample_names <- vapply(fn, get_sample_name, character(1))

message("Merging counts for samples: ", paste(sample_names, collapse = ", "))

count_list <- list()

for (i in seq_along(args)) {
  f <- args[[i]] # get next counts file path
  sname <- sample_names[i]
  
  # Open the counts file as df
  df <- read_tsv(f, comment = "#", show_col_types = FALSE)
  
  if (!"Geneid" %in% colnames(df)) {
    stop(paste("File", f, "does not contain a 'Geneid' column"))
  }
  
  # Selects last column which has counts
  df_sub <- df %>%
    select(Geneid, last_col = ncol(df)) %>%
    rename(!!sname := last_col)
  
  count_list[[sname]] <- df_sub
}

# Combine into a single matrix
merged <- purrr::reduce(count_list, full_join, by = "Geneid") %>%
  arrange(Geneid)

# Replace missing values with zeros
merged[is.na(merged)] <- 0

# Write CSV (matches your declared output filename)
readr::write_csv(merged, outputCountsFile)
