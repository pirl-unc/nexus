tests=(
  test_assembly_long_read_rna
  test_assembly_short_read_rna
)

# Aggregate workflows run several assemblers each, pulling many large images;
# run each in its own pytest invocation and prune Docker between them.
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
