tests=(
  test_alignment_blastp
  test_alignment_bwamem2_1
  test_alignment_bwamem2_2
  test_alignment_diamond_blastp
  test_alignment_minimap2_1
  test_alignment_minimap2_2
  test_alignment_minimap2_dynamic_1
  test_alignment_minimap2_dynamic_2
  test_alignment_star
  test_alignment_ultra
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"
