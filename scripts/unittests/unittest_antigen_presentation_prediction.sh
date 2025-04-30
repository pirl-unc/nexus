pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_antigen_presentation_prediction_mhcflurry2_predict or \
      test_antigen_presentation_prediction_mhcflurry2_predict_scan"
