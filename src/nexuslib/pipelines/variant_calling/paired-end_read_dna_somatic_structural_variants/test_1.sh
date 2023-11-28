#!/bin/bash

# 1. Software paths
NEXTFLOW=/home/jinseok/software/nextflow/nextflow
DELLY2=/home/jinseok/software/delly/delly_v1.1.7_linux_x86_64bit
LUMPY_EXPRESS=/home/jinseok/software/miniconda3/envs/py27/bin/lumpyexpress
LUMPY_EXTRACT_SPLIT_READS_SCRIPTS_FILE=/home/jinseok/software/lumpy/lumpy-sv/scripts/extractSplitReads_BwaMem
BCFTOOLS=/home/jinseok/software/miniconda3/bin/bcftools
PYTHON2=/home/jinseok/software/miniconda3/envs/py27/bin/python
SAMTOOLS=/home/jinseok/software/samtools/samtools-1.16.1/samtools

# 2. Metadata file path
SAMPLES_TSV_FILE=/datastore/lbcfs/collaborations/pirl/workflows/data/raw/metadata/human_illumina_dna_case_control_bam_files.tsv

# 3. File paths
REFERENCE_GENOME_FASTA_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa
LUMPY_CONFIG_FILE=/home/jinseok/software/miniconda3/envs/py27/bin/lumpyexpress.config
DELLY2_EXCLUDE_REGIONS_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/delly2/human.hg38.excl.tsv

# 4. Tool-specific parameters
DELLY2_CALL_PARAMS="--exclude $DELLY2_EXCLUDE_REGIONS_FILE \
                    --map-qual 25 \
                   "
DELLY2_CALL_PARAMS_STR=`echo "$DELLY2_CALL_PARAMS"`

# 5. Nextflow workflow paths
WORK_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/work/variant_calling/paired_end_dna_structural_variants/somatic/
OUTPUT_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/variant_calling/paired_end_dna_structural_variants/somatic/
SCRIPT=paired_end_dna_somatic_structural_variants.nf
CONFIG=../../../nextflow.config

$NEXTFLOW run $SCRIPT -resume \
  -c $CONFIG -profile standard \
  --samples_tsv_file $SAMPLES_TSV_FILE \
  --reference_genome_fasta_file $REFERENCE_GENOME_FASTA_FILE \
  --delly2 $DELLY2 \
  --delly2_call_params "$DELLY2_CALL_PARAMS_STR" \
  --bcftools $BCFTOOLS \
  --python2 $PYTHON2 \
  --lumpy_express $LUMPY_EXPRESS \
  --lumpy_extract_split_reads_script_file $LUMPY_EXTRACT_SPLIT_READS_SCRIPTS_FILE \
  --lumpy_config_file $LUMPY_CONFIG_FILE \
  --samtools $SAMTOOLS \
  --output_dir $OUTPUT_DIR \
  -w $WORK_DIR

