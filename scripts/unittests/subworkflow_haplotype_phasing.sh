tests=(
  test_haplotype_phasing_flair_longshot
  test_haplotype_phasing_hapcut2_whatshap
  test_haplotype_phasing_hiphase
  test_haplotype_phasing_longcallr
  test_haplotype_phasing_longshot
  test_haplotype_phasing_whatshap
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"

