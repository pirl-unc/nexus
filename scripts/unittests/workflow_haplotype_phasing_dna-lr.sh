tests=(
  test_haplotype_phasing_long_read_dna_github_deepvariant
  test_haplotype_phasing_long_read_dna_longshot
  test_haplotype_phasing_long_read_dna_all
  test_haplotype_phasing_long_read_dna_hapcut2_whatshap
  test_haplotype_phasing_long_read_dna_longshot_phaser_only
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"
