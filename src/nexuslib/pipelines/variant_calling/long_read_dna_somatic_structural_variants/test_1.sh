#!/bin/bash

# 1. Software paths
NEXTFLOW=/home/jinseok/software/nextflow/nextflow
SAVANA=/home/jinseok/software/miniconda3/envs/workflows/bin/savana

# 2. Metadata file path
SAMPLES_TSV_FILE=/datastore/lbcfs/collaborations/pirl/workflows/data/raw/metadata/human_pacbio_dna_case_control_bam_files.tsv

# 3. File paths
REFERENCE_GENOME_FASTA_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa

# 4. Tool-specific parameters
SAVANA_PARAMS="--length 30 \
               --mapq 20 \
               --depth 3 \
               "
SAVANA_PARAMS_STR=`echo "$SAVANA_PARAMS"`

# 5. Nextflow workflow paths
WORK_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/work/variant_calling/long_read_dna_somatic_structural_variants/pacbio/
OUTPUT_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/variant_calling/long_read_dna_somatic_structural_variants/pacbio/
SCRIPT=pacbio_dna_somatic_structural_variants.nf
CONFIG=../../../nextflow.config

$NEXTFLOW run $SCRIPT -resume \
  -c $CONFIG -profile standard \
  --samples_tsv_file $SAMPLES_TSV_FILE \
  --reference_genome_fasta_file $REFERENCE_GENOME_FASTA_FILE \
  --savana $SAVANA \
  --savana_params "$SAVANA_PARAMS_STR" \
  --output_dir $OUTPUT_DIR \
  --delete_work_dir false \
  -w $WORK_DIR
