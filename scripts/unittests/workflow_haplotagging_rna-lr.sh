tests=(
  test_haplotagging_long_read_rna_flair_longshot
  test_haplotagging_long_read_rna_longcallr
  test_haplotagging_long_read_rna_all
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"
