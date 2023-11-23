#!/usr/bin/env bash


runLongReadULTRA() {
  TEST_DIR=$1

  # Step 1. Define tools available in environment
  NEXTFLOW=nextflow
  ULTRA=uLTRA
  SAMTOOLS=samtools

  # Step 2. Change directory to current test path
  cd $TEST_DIR

  # Step 3. Define nextflow workflow to test
  SCRIPT=../pipelines/alignment/long_read_rna_ultra/long_read_rna_alignment_ultra.nf

  # Step 4. Create directories
  WORK_DIR=$TEST_DIR/tmp/test_alignment_ultra/work
  INTERMEDIATE=$TEST_DIR/tmp/test_alignment_ultra/intermediate
  OUTPUT_DIR=$TEST_DIR/tmp/test_alignment_ultra/output
  mkdir -p $WORK_DIR
  mkdir -p $INTERMEDIATE
  mkdir -p $OUTPUT_DIR

  # Step 5. Create samples TSV file
  samples_tsv_file=$INTERMEDIATE/test_alignment_ultra_samples.tsv
  printf "%s\t%s\n" "sample_id" "fastq_file" > $samples_tsv_file
  printf "%s\t%s" "sample001" "$TEST_DIR/data/hg38_tp53_variants_rna.fastq.gz" >> $samples_tsv_file

  # Step 6. Define uLTRA parameters
  ULTRA_PARAMS="--isoseq \
               "
  ULTRA_PARAMS_STR=`echo "$ULTRA_PARAMS"`
  REFERENCE_GENOME_FASTA_FILE=$TEST_DIR/data/hg38_chr17_1-8000000.fa
  REFERENCE_GTF_FILE=$TEST_DIR/data/gencode_v41_tp53_annotation.gtf

  # Step 7. Create uLTRA index
  ULTRA_INDEX=$INTERMEDIATE/uLTRA_index/
  mkdir -p $INTERMEDIATE/uLTRA_index/
  $ULTRA index \
    --disable_infer \
    $REFERENCE_GENOME_FASTA_FILE \
    $REFERENCE_GTF_FILE \
    $ULTRA_INDEX

  # Step 7. Run workflow
  $NEXTFLOW run $SCRIPT -resume \
    -c $TEST_DIR/nextflow.config \
    --samples_tsv_file $SAMPLES_TSV_FILE \
    --reference_genome_fasta_file $REFERENCE_GENOME_FASTA_FILE \
    --ultra $ULTRA \
    --ultra_index $ULTRA_INDEX \
    --ultra_params "$ULTRA_PARAMS_STR" \
    --samtools $SAMTOOLS \
    --delete_work_dir false \
    --output_dir $OUTPUT_DIR \
    -w $WORK_DIR
}