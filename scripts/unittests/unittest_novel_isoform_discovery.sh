pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_novel_isoform_discovery_flair or \
      test_novel_isoform_discovery_isoquant or \
      test_novel_isoform_discovery_isoseq"
