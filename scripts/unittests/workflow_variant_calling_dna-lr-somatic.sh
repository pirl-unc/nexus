tests=(
  test_variant_calling_long_read_dna_somatic_github_1
  test_variant_calling_long_read_dna_somatic_2
  test_variant_calling_long_read_dna_somatic_3
)

for test_name in "${tests[@]}"; do
  pytest \
    -s \
    --cov-report=term-missing \
    --cov=nexuslib \
    test/ \
    -k "$test_name"

  docker container prune -f || true
  docker image prune -af || true
  docker builder prune -af || true
done
