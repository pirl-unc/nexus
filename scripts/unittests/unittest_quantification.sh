pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_paired_end_read_rna_quantification_salmon_mapping or \
      test_paired_end_read_rna_quantification_kallisto or \
      test_long_read_rna_quantification_bambu or \
      test_long_read_rna_quantification_kallisto or \
      test_long_read_rna_quantification_liqa or \
      test_long_read_rna_quantification_oarfish"
