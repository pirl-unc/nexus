tests=(
  test_isoform_characterization_espresso
  test_isoform_characterization_flair
  test_isoform_characterization_isoquant
  test_isoform_characterization_isoseq
  test_isoform_characterization_isotools
  test_isoform_characterization_mandalorion
  test_isoform_characterization_rmats
  test_isoform_characterization_sqanti3_fasta
  test_isoform_characterization_sqanti3_gtf
  test_isoform_characterization_talon
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
