pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_haplotype_phasing_long_read_rna_flair_longshot or \
      test_haplotype_phasing_long_read_rna_longcallr or \
      test_haplotype_phasing_long_read_rna_all"
