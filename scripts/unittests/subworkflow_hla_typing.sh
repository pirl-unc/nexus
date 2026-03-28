pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_hla_typing_arcashla or \
      test_hla_typing_hlaprofiler or \
      test_hla_typing_seq2hla"
