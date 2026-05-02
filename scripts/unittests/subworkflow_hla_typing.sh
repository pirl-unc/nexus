tests=(
  test_hla_typing_arcashla
  test_hla_typing_hlaprofiler
  test_hla_typing_seq2hla
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"

