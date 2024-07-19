pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_paired_end_read_rna_hla_typing_arcashla or \
      test_paired_end_read_rna_hla_typing_hlaprofiler or \
      test_docker_paired_end_read_rna_hla_typing_seq2hla"
