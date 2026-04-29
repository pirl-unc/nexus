tests=(
  test_quantification_bambu
  test_quantification_kallisto_lr
  test_quantification_kallisto_pe
  test_quantification_liqa
  test_quantification_oarfish
  test_quantification_salmon_fastq
  test_quantification_transigner
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
