#!/bin/bash

# 1. Software paths
NEXTFLOW=/home/jinseok/software/nextflow/nextflow
GATK4=/home/jinseok/software/gatk4/gatk-4.2.5.0/gatk
PICARD=/home/jinseok/software/picard/picard_2.26.10/picard.jar
STRELKA2=/home/jinseok/software/strelka2/strelka-2.9.10.centos6_x86_64/bin/
SINGULARITY="/usr/bin/singularity"

# 2. Metadata file path
SAMPLES_TSV_FILE=/datastore/lbcfs/collaborations/pirl/workflows/data/raw/metadata/human_illumina_dna_case_control_bam_files.tsv

# 3. File paths
REFERENCE_GENOME_FASTA_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa
GERMLINE_RESOURCE_VCF_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/af-only-gnomad.hg38.vcf
PANEL_OF_NORMALS_VCF_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/1000g_pon.hg38.vcf
KNOWN_VARIANT_SITES_VCF_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/small_exac_common_3.hg38.vcf

# 4. Tool-specific parameters
GATK4_MUTECT2_PARAMS="--germline-resource $GERMLINE_RESOURCE_VCF_FILE \
                      --panel-of-normals $PANEL_OF_NORMALS_VCF_FILE \
                     "
GATK4_MUTECT2_PARAMS_STR=`echo "$GATK4_MUTECT2_PARAMS"`
GATK4_GETPILEUPSUMMARIES_PARAMS="-V $KNOWN_VARIANT_SITES_VCF_FILE \
                                 -L $KNOWN_VARIANT_SITES_VCF_FILE \
                                "
GATK4_GETPILEUPSUMMARIES_PARAMS_STR=`echo "$GATK4_GETPILEUPSUMMARIES_PARAMS"`
STRELKA2_PARAMS=""
STRELKA2_PARAMS_STR=`echo "$STRELKA2_PARAMS"`
DEEPVARIANT_BIN_VERSION="1.4.0"
DEEPVARIANT_BIN_PATH="/opt/deepvariant/bin/run_deepvariant"
DEEPVARIANT_LIB_PATH="/datastore/:/datastore/"
DEEPVARIANT_MODEL_TYPE="WES"
CHROMOSOMES=""
for chr in {1..22}; do
	CHROMOSOMES+=chr"${chr},"
done
CHROMOSOMES+="chrX,"
CHROMOSOMES+="chrY,"
CHROMOSOMES+="chrM"

# 5. Nextflow workflow paths
WORK_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/work/variant_calling/paired_end_dna_small_variants/somatic/human
OUTPUT_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/variant_calling/paired_end_dna_small_variants/somatic/human
SCRIPT=paired_end_human_dna_somatic_small_variants.nf
CONFIG=../../../../nextflow.config

$NEXTFLOW run $SCRIPT -resume \
  -c $CONFIG -profile standard \
  --samples_tsv_file $SAMPLES_TSV_FILE \
  --reference_genome_fasta_file $REFERENCE_GENOME_FASTA_FILE \
  --gatk4 $GATK4 \
  --gatk4_mutect2_params "$GATK4_MUTECT2_PARAMS_STR" \
  --gatk4_getpileupsummaries_params "$GATK4_GETPILEUPSUMMARIES_PARAMS_STR" \
  --picard $PICARD \
  --strelka2 $STRELKA2 \
  --strelka2_params $STRELKA2_PARAMS_STR \
  --singularity $SINGULARITY \
  --deepvariant_bin_path $DEEPVARIANT_BIN_PATH \
  --deepvariant_bin_version $DEEPVARIANT_BIN_VERSION \
  --deepvariant_lib_path $DEEPVARIANT_LIB_PATH \
  --deepvariant_model_type $DEEPVARIANT_MODEL_TYPE \
  --chromosomes $CHROMOSOMES \
  --output_dir $OUTPUT_DIR \
  -w $WORK_DIR
