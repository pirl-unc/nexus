tests=(
  test_hla_typing_arcashla
  test_hla_typing_hlaminer_lr_rna
  test_hla_typing_hlaprofiler
  test_hla_typing_seq2hla
  test_hla_typing_specimmune
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"

