tests=(
  test_utilities_fastq2unalignedbam
  test_utilities_fastqc_1
  test_utilities_fastqc_2
  test_utilities_filter_rnabloom2_transcripts
  test_utilities_pbccs
  test_utilities_sequencing_coverage
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"
