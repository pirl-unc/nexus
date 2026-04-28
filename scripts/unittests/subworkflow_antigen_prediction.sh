tests=(
  test_antigen_prediction_mhcflurry2
  test_antigen_prediction_mhcflurry2_scan
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"
