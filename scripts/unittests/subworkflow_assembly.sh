tests=(
  test_assembly_bookend
  test_assembly_hifiasm
  test_assembly_isonform
  test_assembly_rattle
  test_assembly_rnabloom2
  test_assembly_spades
  test_assembly_stringtie2
  test_assembly_stringtie3
  test_assembly_trinity
)

# Run each tool in its own pytest invocation and prune Docker between them;
# the assemblers pull large images and would otherwise fill the CI disk.
# Note: `test_assembly_rnabloom2` also matches `test_assembly_rnabloom2_clustered`
# (same image), so both run under that entry.
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
