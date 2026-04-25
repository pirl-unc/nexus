pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_haplotyping_flair_longshot or \
      test_haplotyping_longcallr or \
      test_haplotyping_whatshap"
