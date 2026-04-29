tests=(
  test_haplotype_phasing_flair_longshot
  test_haplotype_phasing_hapcut2_whatshap
  test_haplotype_phasing_hiphase
  test_haplotype_phasing_longcallr
  test_haplotype_phasing_longshot
  test_haplotype_phasing_whatshap
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
