tests=(
  test_isoform_characterization_long_read_1
  test_isoform_characterization_long_read_2
  test_isoform_characterization_long_read_3
  test_isoform_characterization_long_read_4
  test_isoform_characterization_short_read
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
