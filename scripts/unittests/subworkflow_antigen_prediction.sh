pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_antigen_prediction_mhcflurry2 or \
      test_antigen_prediction_mhcflurry2_scan"
