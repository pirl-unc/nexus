pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_long_read_dna_variant_calling_cutesv or \
      test_long_read_dna_variant_calling_deepvariant_github or \
      test_long_read_dna_variant_calling_pbsv or \
      test_long_read_dna_variant_calling_savana or \
      test_long_read_dna_variant_calling_severus or \
      test_long_read_dna_variant_calling_sniffles2 or \
      test_long_read_dna_variant_calling_svim or \
      test_long_read_dna_variant_calling_svisionpro"
