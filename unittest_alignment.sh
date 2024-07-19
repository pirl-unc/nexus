pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_long_read_alignment_minimap2_dna and \
      test_long_read_alignment_minimap2_rna and \
      test_long_read_rna_alignment_ultra and \
      test_paired_end_read_dna_alignment_bwa_mem2 and \
      test_paired_end_read_rna_alignment_star"
