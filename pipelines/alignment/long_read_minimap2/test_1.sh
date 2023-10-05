#!/bin/bash

# Tests Oxford Nanopore whole-genome sequencing alignment using minimap2.

# 1. Software paths
NEXTFLOW=/home/jinseok/software/nextflow/nextflow
MINIMAP2=/home/jinseok/software/minimap2/minimap2-2.25/minimap2
SAMTOOLS=/home/jinseok/software/samtools/samtools-1.16.1/samtools

# 2. Metadata file path
SAMPLES_TSV_FILE=/datastore/lbcfs/collaborations/pirl/workflows/data/raw/metadata/human_ont_dna_fastq_files.tsv

# 3. File paths
REFERENCE_GENOME_FASTA_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa

# 4. Tool-specific parameters
MINIMAP2_PARAMS="-ax map-ont \
                 --cs \
                 --eqx \
                 -Y \
                 -L"
MINIMAP2_PARAMS_STR=`echo "$MINIMAP2_PARAMS"`

# 5. Nextflow workflow paths
WORK_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/work/alignment/long_read_minimap2/ont_dna_minimap2/
OUTPUT_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/alignment/long_read_minimap2/ont_dna_minimap2/
SCRIPT=long_read_alignment_minimap2.nf
CONFIG=../../nextflow.config

$NEXTFLOW run $SCRIPT -resume \
  -c $CONFIG -profile standard \
  --samples_tsv_file $SAMPLES_TSV_FILE \
  --reference_genome_fasta_file $REFERENCE_GENOME_FASTA_FILE \
  --minimap2 $MINIMAP2 \
  --minimap2_params "$MINIMAP2_PARAMS_STR" \
  --samtools $SAMTOOLS \
  --platform_tag "ont" \
  --delete_work_dir true \
  --output_dir $OUTPUT_DIR \
  -w $WORK_DIR

