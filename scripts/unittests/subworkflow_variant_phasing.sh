pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_variant_phasing_hiphase or \
      test_variant_phasing_whatshap"
