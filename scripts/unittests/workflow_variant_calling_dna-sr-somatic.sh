pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_variant_calling_short_read_dna_somatic_github_1 or \
      test_variant_calling_short_read_dna_somatic_2 or \
      test_variant_calling_short_read_dna_somatic_3 or \
      test_variant_calling_short_read_dna_somatic_4"
