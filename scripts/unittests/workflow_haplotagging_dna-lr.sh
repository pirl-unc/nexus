tests=(
  test_haplotagging_long_read_dna_github_deepvariant
  test_haplotagging_long_read_dna_longshot
  test_haplotagging_long_read_dna_all
  test_haplotagging_long_read_dna_hapcut2_whatshap
  test_haplotagging_long_read_dna_longphase
  test_haplotagging_long_read_dna_longshot_phaser_only
)

for test_name in "${tests[@]}"; do
  pytest \
    -s \
    --cov-report=term-missing \
    --cov=nexuslib \
    test/ \
    -k "$test_name"
done
