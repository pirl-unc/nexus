pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "not test_long_read_dna_small_variants_docker and \
      not test_long_read_dna_structural_variants_local and \
      not test_paired_end_read_dna_somatic_small_variants_human_docker"