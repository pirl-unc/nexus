pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_haplotype_phasing_flair_longshot or \
      test_haplotype_phasing_hapcut2_whatshap or \
      test_haplotype_phasing_hiphase or \
      test_haplotype_phasing_longcallr or \
      test_haplotype_phasing_longshot or \
      test_haplotype_phasing_whatshap"
