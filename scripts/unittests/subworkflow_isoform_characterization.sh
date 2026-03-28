pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_isoform_characterization_espresso or \
      test_isoform_characterization_flair or \
      test_isoform_characterization_isoquant or \
      test_isoform_characterization_isoseq or \
      test_isoform_characterization_isotools or \
      test_isoform_characterization_mandalorion or \
      test_isoform_characterization_rmats or \
      test_isoform_characterization_sqanti3_fasta or \
      test_isoform_characterization_sqanti3_gtf or \
      test_isoform_characterization_talon"
