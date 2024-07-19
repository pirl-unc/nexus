pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_long_read_dna_variant_calling_cutesv or \
      test_long_read_dna_variant_calling_deepvariant or \
      test_long_read_dna_variant_calling_pbsv or \
      test_long_read_dna_variant_calling_savana or \
      test_long_read_dna_variant_calling_severus"
