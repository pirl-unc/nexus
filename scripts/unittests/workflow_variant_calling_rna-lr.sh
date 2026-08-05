tests=(
  test_variant_calling_long_read_rna_1
  test_variant_calling_long_read_rna_github_2
)

# Aggregate RNA long-read variant-calling workflow runs several callers
# (incl. DeepVariant via de-souza); run each in its own pytest invocation and
# prune Docker between them. The `*_local_*` variants are excluded (they need
# local-only data not present on CI).
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
