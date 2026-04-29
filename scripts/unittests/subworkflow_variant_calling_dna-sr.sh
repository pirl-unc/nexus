tests=(
  test_variant_calling_clair3_2
  test_variant_calling_clairs_2
  test_variant_calling_deepsomatic_github_2
  test_variant_calling_deepvariant_github_2
  test_variant_calling_delly2_sr_germline
  test_variant_calling_delly2_sr_somatic
  test_variant_calling_dysgu_germline_pe
  test_variant_calling_dysgu_somatic_pe
  test_variant_calling_gridss2_germline
  test_variant_calling_gridss2_somatic
  test_variant_calling_haplotypecaller
  test_variant_calling_lumpy_germline
  test_variant_calling_lumpy_somatic
  test_variant_calling_manta_germline
  test_variant_calling_manta_somatic
  test_variant_calling_mutect2_1
  test_variant_calling_mutect2_2
  test_variant_calling_octopus_germline
  test_variant_calling_octopus_somatic
  test_variant_calling_pindel
  test_variant_calling_sequenza
  test_variant_calling_strelka2_germline
  test_variant_calling_strelka2_somatic
  test_variant_calling_svaba
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
