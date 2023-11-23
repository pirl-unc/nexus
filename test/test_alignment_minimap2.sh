#!/usr/bin/env bash


runLongReadMinimap2() {
  TEST_DIR=$1

  # Step 1. Define tools available in environment
  NEXTFLOW=nextflow
  MINIMAP2=minimap2
  SAMTOOLS=samtools

  # Step 2. Change directory to current test path
  cd $TEST_DIR

  # Step 3. Define nextflow workflow to test
  SCRIPT=../pipelines/alignment/long_read_minimap2/long_read_alignment_minimap2.nf

  # Step 4. Create directories
  WORK_DIR=$TEST_DIR/tmp/test_alignment_minimap2/work
  INTERMEDIATE=$TEST_DIR/tmp/test_alignment_minimap2/intermediate
  OUTPUT_DIR=$TEST_DIR/tmp/test_alignment_minimap2/output
  mkdir -p $WORK_DIR
  mkdir -p $INTERMEDIATE
  mkdir -p $OUTPUT_DIR

  # Step 5. Create samples TSV file
  samples_tsv_file=$INTERMEDIATE/test_alignment_minimap2_samples.tsv
  printf "%s\t%s\n" "sample_id" "fastq_file" > $samples_tsv_file
  printf "%s\t%s" "sample001" "$TEST_DIR/data/hg38_tp53_variants_dna.fastq.gz" >> $samples_tsv_file

  # Step 6. Define minimap2 parameters
  MINIMAP2_PARAMS="-ax map-hifi --cs --eqx -Y -L"
  MINIMAP2_PARAMS_STR=`echo "$MINIMAP2_PARAMS"`

  # Step 7. Run workflow
  nextflow run $SCRIPT -resume \
    -c $TEST_DIR/nextflow.config \
    --samples_tsv_file $samples_tsv_file \
    --reference_genome_fasta_file $TEST_DIR/data/hg38_chr17_1-8000000.fa \
    --minimap2 $MINIMAP2 \
    --minimap2_params "$MINIMAP2_PARAMS_STR" \
    --samtools $SAMTOOLS \
    --platform_tag "pacbio" \
    --delete_work_dir false \
    --output_dir $OUTPUT_DIR \
    -w $WORK_DIR
}