pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_haplotype_phasing_long_read_dna_github_deepvariant or \
      test_haplotype_phasing_long_read_dna_longshot or \
      test_haplotype_phasing_long_read_dna_all or \
      test_haplotype_phasing_long_read_dna_hapcut2_whatshap or \
      test_haplotype_phasing_long_read_dna_longshot_phaser_only"
