pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_transcriptome_assembly_rnabloom2 or \
      test_transcriptome_assembly_stringtie2"
