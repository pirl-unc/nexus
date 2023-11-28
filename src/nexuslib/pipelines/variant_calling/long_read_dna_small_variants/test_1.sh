#!/bin/bash

# 1. Software paths
NEXTFLOW=/home/jinseok/software/nextflow/nextflow
SINGULARITY="/usr/bin/singularity"

# 2. Metadata file path
SAMPLES_TSV_FILE=/datastore/lbcfs/collaborations/pirl/workflows/data/raw/metadata/human_pacbio_dna_bam_files.tsv

# 3. File paths
REFERENCE_GENOME_FASTA_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa

# 4. Tool-specific parameters
DEEPVARIANT_BIN_VERSION="1.4.0"
DEEPVARIANT_BIN_PATH="/opt/deepvariant/bin/run_deepvariant"
DEEPVARIANT_LIB_PATH="/datastore/:/datastore/"
DEEPVARIANT_MODEL_TYPE="PACBIO"

# 5. Nextflow workflow paths
WORK_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/work/variant_calling/long_read_dna_small_variants/pacbio/
OUTPUT_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/variant_calling/long_read_dna_small_variants/pacbio/
SCRIPT=pacbio_dna_small_variants.nf
CONFIG=../../../nextflow.config

$NEXTFLOW run $SCRIPT -resume \
  -c $CONFIG -profile standard \
  --samples_tsv_file $SAMPLES_TSV_FILE \
  --reference_genome_fasta_file $REFERENCE_GENOME_FASTA_FILE \
  --singularity $SINGULARITY \
  --deepvariant_bin_path $DEEPVARIANT_BIN_PATH \
  --deepvariant_bin_version $DEEPVARIANT_BIN_VERSION \
  --deepvariant_lib_path $DEEPVARIANT_LIB_PATH \
  --deepvariant_model_type $DEEPVARIANT_MODEL_TYPE \
  --output_dir $OUTPUT_DIR \
  --delete_work_dir false \
  -w $WORK_DIR
