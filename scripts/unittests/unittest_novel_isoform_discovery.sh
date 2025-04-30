pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_novel_isoform_discovery_espresso or \
      test_novel_isoform_discovery_flair or \
      test_novel_isoform_discovery_isoquant or \
      test_novel_isoform_discovery_isoseq or \
      test_novel_isoform_discovery_isotools or \
      test_novel_isoform_discovery_mandalorion or \
      test_novel_isoform_discovery_sqanti3_fastq_mode or \
      test_novel_isoform_discovery_sqanti3_gtf_mode_flair or \
      test_novel_isoform_discovery_sqanti3_gtf_mode_isoquant or \
      test_novel_isoform_discovery_talon"
