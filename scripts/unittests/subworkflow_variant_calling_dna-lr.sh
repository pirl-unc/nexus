tests=(
  test_variant_calling_clair3_1
  test_variant_calling_clairs_1
  test_variant_calling_cutesv
  test_variant_calling_deepsomatic_github_1
  test_variant_calling_deepvariant_github_1
  test_variant_calling_delly2_lr_germline
  test_variant_calling_delly2_lr_somatic
  test_variant_calling_dysgu_somatic_lr
  test_variant_calling_hificnv
  test_variant_calling_longshot
  test_variant_calling_nanocaller_1
  test_variant_calling_nanocaller_2
  test_variant_calling_nanomonsv
  test_variant_calling_nanovar_1
  test_variant_calling_nanovar_2
  test_variant_calling_pbsv
  test_variant_calling_savana
  test_variant_calling_severus
  test_variant_calling_sniffles2
  test_variant_calling_svisionpro
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"
