pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_utilities_fastq2unalignedbam or \
      test_utilities_fastqc_1 or \
      test_utilities_fastqc_2 or \
      test_utilities_filter_rnabloom2_transcripts or \
      test_utilities_pbccs or \
      test_utilities_sequencing_coverage"
