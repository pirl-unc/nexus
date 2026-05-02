tests=(
  test_variant_calling_clair3rna
  test_variant_calling_de_souza_github
  test_variant_calling_longgf
  test_variant_calling_pbfusion
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"
