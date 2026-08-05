tests=(
  test_variant_calling_clair3rna
  test_variant_calling_de_souza_github
  test_variant_calling_isolaser
  test_variant_calling_longgf
  test_variant_calling_pbfusion
)

# Run each tool in its own pytest invocation and prune Docker between them;
# de-souza (DeepVariant) and isolaser (GATK) pull large images that would
# otherwise fill the CI disk if run together.
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
