#!/bin/bash

# 1. Software paths
NEXTFLOW=/home/jinseok/software/nextflow/nextflow
SALMON=/home/jinseok/software/salmon/salmon-1.8.0_linux_x86_64/bin/salmon

# 2. Metadata file path
SAMPLES_TSV_FILE=/datastore/lbcfs/collaborations/pirl/workflows/data/raw/metadata/human_illumina_rna_fastq_files.tsv

# 3. File paths
GTF_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf

# 4. Tool-specific parameters
SALMON_INDEX_PARAMS="--gencode \
                    "
SALMON_INDEX_PARAMS_STR=`echo "$SALMON_INDEX_PARAMS"`
SALMON_QUANT_PARAMS="--libType IU \
                     --geneMap $GTF_FILE \
                     --seqBias \
                     --gcBias \
                     --posBias \
                    "
SALMON_QUANT_PARAMS_STR=`echo "$SALMON_QUANT_PARAMS"`

# 5. Nextflow workflow paths
WORK_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/work/quantification/paired_end_read_rna/salmon_mapping_mode/
OUTPUT_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/quantification/paired_end_read_rna/salmon_mapping_mode/
SCRIPT=paired_end_read_rna_quantification_salmon_mapping_mode.nf
CONFIG=../../../nextflow.config

$NEXTFLOW run $SCRIPT -resume \
  -c $CONFIG -profile standard \
  --samples_tsv_file $SAMPLES_TSV_FILE \
  --salmon $SALMON \
  --salmon_index_params "$SALMON_INDEX_PARAMS_STR" \
  --salmon_quant_params "$SALMON_QUANT_PARAMS_STR" \
  --output_dir $OUTPUT_DIR \
  --delete_work_dir false \
  -w $WORK_DIR
