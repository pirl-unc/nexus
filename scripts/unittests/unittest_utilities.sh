pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_fastq_to_unaligned_bam or \
      test_fastqc_single_end_reads or \
      test_fastqc_paired_end_reads or \
      test_sequencing_coverage"
