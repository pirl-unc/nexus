#!/bin/bash

# Tests PacBio RNA sequencing alignment using ultra.

# 1. Software paths
NEXTFLOW=/home/jinseok/software/nextflow/nextflow
ULTRA=/home/jinseok/software/miniconda3/bin/uLTRA
SAMTOOLS=/home/jinseok/software/samtools/samtools-1.16.1/samtools

# 2. Metadata file path
SAMPLES_TSV_FILE=/datastore/lbcfs/collaborations/pirl/workflows/data/raw/metadata/human_pacbio_rna_fastq_files.tsv

# 3. File paths
REFERENCE_GENOME_FASTA_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa
ULTRA_INDEX=/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/ultra/hg38_index/

# 4. Tool-specific parameters
ULTRA_PARAMS="--isoseq \
             "
ULTRA_PARAMS_STR=`echo "$ULTRA_PARAMS"`

# 5. Nextflow workflow paths
WORK_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/work/alignment/long_read_rna_ultra/pacbio_rna_ultra
OUTPUT_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/alignment/long_read_rna_ultra/pacbio_rna_ultra
SCRIPT=long_read_rna_alignment_ultra.nf
CONFIG=../../nextflow.config

$NEXTFLOW run $SCRIPT -resume \
  -c $CONFIG -profile standard \
  --samples_tsv_file $SAMPLES_TSV_FILE \
  --reference_genome_fasta_file $REFERENCE_GENOME_FASTA_FILE \
  --ultra $ULTRA \
  --ultra_index $ULTRA_INDEX \
  --ultra_params "$ULTRA_PARAMS_STR" \
  --samtools $SAMTOOLS \
  --delete_work_dir false \
  --output_dir $OUTPUT_DIR \
  -w $WORK_DIR
