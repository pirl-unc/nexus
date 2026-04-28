tests=(
  test_variant_calling_long_read_dna_germline_1
  test_variant_calling_long_read_dna_germline_2
  test_variant_calling_long_read_dna_germline_github_3
  test_variant_calling_long_read_dna_somatic_github_1
  test_variant_calling_long_read_dna_somatic_2
  test_variant_calling_long_read_dna_somatic_3
  test_variant_calling_short_read_dna_germline_1
  test_variant_calling_short_read_dna_germline_2
  test_variant_calling_short_read_dna_germline_github_3
  test_variant_calling_short_read_dna_somatic_github_1
  test_variant_calling_short_read_dna_somatic_2
  test_variant_calling_short_read_dna_somatic_3
  test_variant_calling_short_read_dna_somatic_4
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"
