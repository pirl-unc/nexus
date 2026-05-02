tests=(
  test_haplotagging_short_read_dna_github_deepvariant
  test_haplotagging_short_read_dna_haplotypecaller
  test_haplotagging_short_read_dna_clair3
  test_haplotagging_short_read_dna_strelka2_germline
  test_haplotagging_short_read_dna_all_github
)

for test_name in "${tests[@]}"; do
  pytest \
    -s \
    --cov-report=term-missing \
    --cov=nexuslib \
    test/ \
    -k "$test_name"
done
