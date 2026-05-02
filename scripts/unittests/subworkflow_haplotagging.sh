tests=(
  test_haplotagging_flair_longshot
  test_haplotagging_hapcut2_whatshap
  test_haplotagging_hiphase
  test_haplotagging_longcallr
  test_haplotagging_longphase
  test_haplotagging_longshot
  test_haplotagging_margin
  test_haplotagging_whatshap
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
