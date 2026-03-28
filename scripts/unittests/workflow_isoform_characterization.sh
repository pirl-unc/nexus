pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_isoform_characterization_long_read_1 or \
      test_isoform_characterization_long_read_2 or \
      test_isoform_characterization_long_read_3 or \
      test_isoform_characterization_short_read"
