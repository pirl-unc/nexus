pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_assembly_hifiasm or \
      test_assembly_rnabloom2 or \
      test_assembly_stringtie2"
