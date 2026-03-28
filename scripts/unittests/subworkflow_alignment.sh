pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_alignment_blastp or \
      test_alignment_bwamem2_1 or \
      test_alignment_bwamem2_2 or \
      test_alignment_diamond_blastp or \
      test_alignment_minimap2_1 or \
      test_alignment_minimap2_2 or \
      test_alignment_minimap2_dynamic_1 or \
      test_alignment_minimap2_dynamic_2 or \
      test_alignment_star or \
      test_alignment_ultra"
