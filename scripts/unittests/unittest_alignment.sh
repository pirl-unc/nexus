pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_long_read_alignment_minimap2_dna or \
      test_long_read_alignment_minimap2_rna or \
      test_long_read_rna_alignment_ultra or \
      test_paired_end_read_dna_alignment_bwa_mem2 or \
      test_paired_end_read_rna_alignment_star"
