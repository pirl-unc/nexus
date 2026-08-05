tests=(
  test_alignment_blastp
  test_alignment_bwamem2_1
  test_alignment_bwamem2_2
  test_alignment_bwamem2_3
  test_alignment_diamond_blastp
  test_alignment_minimap2_1
  test_alignment_minimap2_2
  test_alignment_minimap2_3
  test_alignment_minimap2_dynamic_1
  test_alignment_minimap2_dynamic_2
  test_alignment_star
  test_alignment_ultra
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
