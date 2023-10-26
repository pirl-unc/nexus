#!/bin/bash

# 1. Software paths
NEXTFLOW=/home/jinseok/software/nextflow/nextflow
SNIFFLES2=/home/jinseok/software/miniconda3/envs/exacto/bin/sniffles
SVIM=/home/jinseok/software/miniconda3/envs/exacto/bin/svim
CUTESV=/home/jinseok/software/miniconda3/envs/exacto/bin/cuteSV
PBSV=/home/jinseok/software/miniconda3/envs/exacto/bin/pbsv

# 2. Metadata file path
SAMPLES_TSV_FILE=/datastore/lbcfs/collaborations/pirl/workflows/data/raw/metadata/human_pacbio_dna_bam_files.tsv

# 3. File paths
REFERENCE_GENOME_FASTA_FILE=/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa

# 4. Tool-specific parameters
SNIFFLES2_PARAMS="--minsupport 3 \
                  --minsvlen 30 \
                  --mapq 20 \
                  --output-rnames \
                 "
SNIFFLES2_PARAMS_STR=`echo "$SNIFFLES2_PARAMS"`
PBSV_DISCOVER_PARAMS="--ccs \
                      --min-gap-comp-id-perc 97 \
                      --min-mapq 20 \
                     "
PBSV_DISCOVER_PARAMS_STR=`echo "$PBSV_DISCOVER_PARAMS"`
PBSV_CALL_PARAMS="--ccs \
                  --call-min-reads-per-strand-all-samples 0 \
                  --call-min-read-perc-one-sample 10 \
                  --call-min-reads-all-samples 3 \
                  --call-min-reads-one-sample 3 \
                 "
PBSV_CALL_PARAMS_STR=`echo "$PBSV_CALL_PARAMS"`
SVIM_PARAMS="--min_mapq 20 \
             --min_sv_size 30 \
             --insertion_sequences \
             --read_names \
             --zmws \
            "
SVIM_PARAMS_STR=`echo "$SVIM_PARAMS"`
CUTESV_PARAMS="--max_cluster_bias_INS 1000 \
               --diff_ratio_merging_INS	0.9 \
               --max_cluster_bias_DEL	1000 \
               --diff_ratio_merging_DEL	0.5 \
               --min_support 3 \
               --min_mapq	20 \
               --min_size	30 \
               --max_size -1 \
               --report_readid \
               --genotype \
              "
CUTESV_PARAMS_STR=`echo "$CUTESV_PARAMS"`
EXACTO_CONVERT_PARAMS="--sequencing-platform pacbio \
                      "
EXACTO_CONVERT_PARAMS_STR=`echo "$EXACTO_CONVERT_PARAMS"`

# 5. Nextflow workflow paths
WORK_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/work/variant_calling/long_read_dna_structural_variants/pacbio/
OUTPUT_DIR=/datastore/lbcfs/collaborations/pirl/workflows/data/processed/variant_calling/long_read_dna_structural_variants/pacbio/
SCRIPT=pacbio_dna_structural_variants.nf
CONFIG=../../../nextflow.config

$NEXTFLOW run $SCRIPT -resume \
  -c $CONFIG -profile standard \
  --samples_tsv_file $SAMPLES_TSV_FILE \
  --reference_genome_fasta_file $REFERENCE_GENOME_FASTA_FILE \
  --sniffles2 $SNIFFLES2 \
  --sniffles2_params "$SNIFFLES2_PARAMS_STR" \
  --pbsv $PBSV \
  --pbsv_discover_params "$PBSV_DISCOVER_PARAMS_STR" \
  --pbsv_call_params "$PBSV_CALL_PARAMS_STR" \
  --svim $SVIM \
  --svim_params "$SVIM_PARAMS_STR" \
  --cutesv $CUTESV \
  --cutesv_params "$CUTESV_PARAMS_STR" \
  --exacto $EXACTO \
  --exacto_convert_params "$EXACTO_CONVERT_PARAMS_STR" \
  --output_dir $OUTPUT_DIR \
  --delete_work_dir false \
  -w $WORK_DIR
