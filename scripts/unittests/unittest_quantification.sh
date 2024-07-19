pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_paired_end_read_rna_quantification_salmon_mapping"
