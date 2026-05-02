tests=(
  test_isoform_characterization_long_read_1
  test_isoform_characterization_long_read_2
  test_isoform_characterization_long_read_3
  test_isoform_characterization_short_read
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"
