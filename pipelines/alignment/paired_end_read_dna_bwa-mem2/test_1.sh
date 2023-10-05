#!/bin/bash

# Tests Illumina paired-end DNA sequencing alignment using bwa-mem2.

# 1. Software paths
NEXTFLOW=/home/jinseok/software/nextflow/nextflow
BWA_MEM2=/home/jinseok/software/bwa_mem2/bwa-mem2-2.2.1_x64-linux/bwa-mem2
SAMTOOLS=/home/jinseok/software/samtools/samtools-1.16.1/samtools
ABRA2=/home/jinseok/software/abra2/abra2-2.23.jar
GATK4=/home/jinseok/software/gatk4/gatk-4.2.5.0/gatk

# 2. Metadata file path
SAMPLES_TSV_FILE=/datastore/lbcfs/collaborations/pirl/workflows/data/raw/metadata/human_illumina_dna_fastq_files.tsv

# 3. File paths
REFERENCE_GENOME_FASTA_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa
ABRA2_TARGETS_BED_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/abra2/gencode-v41-annotation-abra2-exon-targets.bed
DBSNP_SNP_VCF_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/dbsnp_146.hg38.vcf
THOUSAND_GENOMES_VCF_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/1000G_phase1.snps.high_confidence.hg38.vcf
MILLS_AND_1000G_GOLD_STANDARD_VCF_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/Mills_and_1000G_gold_standard.indels.hg38.vcf
KNOWN_INDELS_VCF_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/Homo_sapiens_assembly38.known_indels.vcf

# 4. Tool-specific parameters
GATK4_BASERECALIBRATOR_PARAMS="--known-sites $DBSNP_SNP_VCF_FILE \
                               --known-sites $THOUSAND_GENOMES_VCF_FILE \
                               --known-sites $MILLS_AND_1000G_GOLD_STANDARD_VCF_FILE \
                               --known-sites $KNOWN_INDELS_VCF_FILE \
                              "
GATK4_BASERECALIBRATOR_PARAMS_STR=`echo "$GATK4_BASERECALIBRATOR_PARAMS"`
CHROMOSOMES=""
for chr in {1..22}; do
	CHROMOSOMES+=chr"${chr},"
done
CHROMOSOMES+="chrX,"
CHROMOSOMES+="chrY,"
CHROMOSOMES+="chrM"

# 5. Nextflow workflow paths
WORK_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/work/alignment/paired_end_read_dna_bwa-mem2/
OUTPUT_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/alignment/paired_end_read_dna_bwa-mem2/
SCRIPT=paired_end_read_dna_alignment_bwa-mem2.nf
CONFIG=../../nextflow.config

$NEXTFLOW run $SCRIPT -resume \
  -c $CONFIG -profile standard \
  --samples_tsv_file $SAMPLES_TSV_FILE \
  --reference_genome_fasta_file $REFERENCE_GENOME_FASTA_FILE \
  --bwa_mem2 $BWA_MEM2 \
  --abra2 $ABRA2 \
  --gatk4 $GATK4 \
  --gatk4_baserecalibrator_params "$GATK4_BASERECALIBRATOR_PARAMS_STR" \
  --samtools $SAMTOOLS \
  --abra2_targets_bed_file $ABRA2_TARGETS_BED_FILE \
  --chromosomes "$CHROMOSOMES" \
  --output_dir $OUTPUT_DIR \
  --delete_work_dir true \
  -w $WORK_DIR
