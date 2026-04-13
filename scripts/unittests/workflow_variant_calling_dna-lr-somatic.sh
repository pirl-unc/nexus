pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_variant_calling_long_read_dna_somatic_github_1 or \
      test_variant_calling_long_read_dna_somatic_2 or \
      test_variant_calling_long_read_dna_somatic_3"
