pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_variant_calling_long_read_germline_1 or \
      test_variant_calling_long_read_germline_2 or \
      test_variant_calling_long_read_germline_github_3 or \
      test_variant_calling_long_read_somatic_github_1 or \
      test_variant_calling_long_read_somatic_2 or \
      test_variant_calling_long_read_somatic_3 or \
      test_variant_calling_short_read_germline_1 or \
      test_variant_calling_short_read_germline_2 or \
      test_variant_calling_short_read_germline_github_3 or \
      test_variant_calling_short_read_somatic_github_1 or \
      test_variant_calling_short_read_somatic_2 or \
      test_variant_calling_short_read_somatic_3 or \
      test_variant_calling_short_read_somatic_4"
