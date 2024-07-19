pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_long_read_dna_variant_phasing_hiphase or \
      test_long_read_dna_variant_phasing_whatshap"
