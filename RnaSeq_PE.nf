#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Users can override default settings using the nextflow.config file

// Input files
params.reads    = 		"${baseDir}/out/clean/*_{R1,R2}.clean.fastq.gz"
params.metadata = 		"${baseDir}/metadata.csv"

// Script locations
params.deseq2script = 	"${baseDir}/scripts/DESeq2.R"
params.ctMatrixScript = "${baseDir}/scripts/makeCountsMatrix.R"

// Genome files
params.starGenomeDir = 	"${baseDir}/genomes/star_idx_chrM"
params.fc_gtf = 		"${baseDir}/genomes/featureCounts/gc48.chrM.gtf"

// Output directory organization
params.starOutDir = 	"${baseDir}/analysis/star_output"
params.starBamDir = 	"${baseDir}/analysis/star_bam"
params.umiOutDir = 		"${baseDir}/analysis/umi_dedup"
params.umiCountOutDir = "${baseDir}/analysis/umi_count"
params.countsOutDir = 	"${baseDir}/analysis/counts"
params.deOutDir = 		"${baseDir}/analysis/de"
params.fcOutDir = 		"${baseDir}/analysis/featureCounts"

// featureCounts parameters
params.fc_feature = 'exon'
params.fc_attr = 'gene_id'
params.fc_strand = '0'
params.fc_extra = ''

// Count Matrix File Name (must be .csv)
params.countsMatrix = "countsMatrix.csv"

// UMI deduplicate params
params.umi_tag = 'RX'
params.paired_end = '--paired'
params.umi_spliceduniq = '--spliced-is-unique'
params.umi_dedup_opts = ''
params.use_umi_dedup = true

//*Fastq file channel *//
Channel
    .fromPath("$params.metadata")
    .splitCsv(header:true)
    .map { row ->
		def fq1 = file(row.fastq_1)
        def fq2 = file(row.fastq_2)
		def smp = fq1.getName().tokenize('_')[0]
		

        if (!fq1.exists()) {
            error "Missing: $fq1"
        } else {
            println "File found: $fq1"
        }

        if (!fq2.exists()) {
            error "Missing: $fq2"
        } else {
            println "File found: $fq2"
        }
		[smp, fq1, fq2]
    }
    .set { SAMPLES }

//*Get differential expression analysis groups*//
Channel
    .fromPath("$params.metadata")
    .splitCsv(header:true)
    .map { row ->
		def group = file(row.condition)
        def smp = file(row.fastq_1).getName().tokenize('_')[0]
		[smp, group]
    }
	.set { DE_GROUP }

//* Read alignment - splicing aware *//
process ALIGN_WITH_STAR {
    tag { sample_id }
    
	input:
    tuple val(sample_id), path(fq1), path(fq2)

	publishDir "$params.starOutDir", pattern: '*.Log.final.out', mode: 'copy'
	publishDir "$params.starOutDir", pattern: '*.Log.out', mode: 'copy'
	publishDir "$params.starOutDir", pattern: '*.Log.progress.out', mode: 'copy'
	publishDir "$params.starOutDir", pattern: '*.SJ.progress.tab', mode: 'copy'
	publishDir "$params.starBamDir", pattern: '*.bam', mode: 'copy'
	
    output:
    tuple val(sample_id), path("*.Aligned.out.bam"), emit: sampleid_and_genome_bam_path
	tuple val(sample_id), path("*.Aligned.toTranscriptome.out.bam"), emit: sampleid_and_trx_bam_path
	path ("*"), emit: star_out

    script:
    """
    STAR \
        --runThreadN 4 \
        --genomeDir ${params.starGenomeDir} \
        --readFilesIn ${fq1} ${fq2} \
        --readFilesCommand zcat \
		--outSAMtype BAM Unsorted \
        --outFileNamePrefix ${sample_id}. \
		--quantMode TranscriptomeSAM
    """    
}

process FEATURECOUNTS {
  tag "$sample"
  label 'counts_small'
  publishDir "${params.fcOutDir}", mode: 'copy'

  cpus 4
  memory '8 GB'
  errorStrategy 'terminate'

  input:
  tuple val(sample), path(bam)   // genome-aligned BAM (unsorted ok)
  path genome_gtf

  output:
  tuple val(sample), path("${sample}.featureCounts.bam"),     emit: FC_BAM
  path("${sample}.featureCounts.bam.bai"),                    emit: FC_BAI
  path("${sample}.featureCounts.txt"),                        emit: FC_TABLE
  path("${sample}.featureCounts.summary"),                    emit: FC_SUMMARY
  // optional: keep coord-sorted input for downstream
  path("${sample}.coord.bam"),                                optional: true
  path("${sample}.coord.bam.bai"),                            optional: true

  script:
  """
  set -euo pipefail

  test -s "${genome_gtf}" || { echo "[ERROR] GTF not found: ${genome_gtf}" >&2; exit 1; }
  test -s "${bam}" || { echo "[ERROR] Input BAM missing: ${bam}" >&2; exit 1; }

  # Ensure coordinate sort (safe if already sorted)
  samtools sort -@ ${task.cpus} -o ${sample}.coord.bam ${bam}
  samtools index -@ ${task.cpus} ${sample}.coord.bam

  # Tag reads with gene assignment (XT/XS) and output counts
  featureCounts \\
    -p \\
	-T ${task.cpus} \\
    -a ${genome_gtf} \\
    -t ${params.fc_feature} \\
	-g ${params.fc_attr} \\
	-s ${params.fc_strand} \\
	-R BAM \\
    -o ${sample}.featureCounts.txt \\
    ${params.fc_extra} \\
    ${sample}.coord.bam

  # featureCounts writes '<input>.featureCounts.bam' next to input
  mv ${sample}.coord.bam.featureCounts.bam ${sample}.featureCounts.bam
  samtools index -@ ${task.cpus} ${sample}.featureCounts.bam

  # Normalize summary name if present
  if [ -f "${sample}.featureCounts.txt.summary" ]; then
    mv "${sample}.featureCounts.txt.summary" "${sample}.featureCounts.summary"
  fi
  """
}

//* From counts file list - make a single counts matrix.
// Counts files must have gene colum header "Geneid" *//
process MAKE_COUNT_MATRIX {

  publishDir "${params.countsOutDir}", mode: 'copy'

  input:
  path(featureCountFileList)
  path(makeCountMatrixRscript)

  output:
  path("*.csv"), emit: MATRIX

  script:
  """
  set -euo pipefail
  Rscript ${makeCountMatrixRscript} ${featureCountFileList} ${params.countsMatrix}
  """
}

//*Run DESeq2 on the counts matrix *//
process RUN_DESEQ2 {
	tag "Run DESeq2"
	stageInMode 'copy' 

	publishDir "${params.deOutDir}", mode: "copy"

	input:
	path(deseq2script)
	path(countsMatrix)
	path(metadata)

	output:
	path("*")

	script:
	"""
	Rscript ${deseq2script} ${countsMatrix} ${metadata} 
	"""
}

workflow {
    // --- Step 1: Run STAR mapping ---
	STAR_OUT = ALIGN_WITH_STAR( SAMPLES )

	// --- Step 2: Run FeatureCounts ---
	FC_OUT = FEATURECOUNTS(STAR_OUT.sampleid_and_genome_bam_path, params.fc_gtf)

	// --- Step 3: Make a counts matrix for DESeq2 ---
	COUNT_MATRIX = MAKE_COUNT_MATRIX(FC_OUT.FC_TABLE.collect(), params.ctMatrixScript)

	// --- Step 4: Run DEseq2 ---
	RUN_DESEQ2(params.deseq2script, COUNT_MATRIX.MATRIX, params.metadata)
}