pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_quantification_long_read_1 or \
      test_quantification_long_read_2 or \
      test_quantification_short_read"
