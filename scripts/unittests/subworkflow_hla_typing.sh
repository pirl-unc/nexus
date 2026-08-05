tests=(
  test_hla_typing_arcashla
  test_hla_typing_fufihla
  test_hla_typing_hlaminer_lr_dna
  test_hla_typing_hlaminer_lr_rna
  test_hla_typing_hlaminer_sr_dna
  test_hla_typing_hlaminer_sr_rna
  test_hla_typing_hlaprofiler
  test_hla_typing_optitype
  test_hla_typing_seq2hla
  test_hla_typing_spechla
  test_hla_typing_spechla_lr_dna
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

