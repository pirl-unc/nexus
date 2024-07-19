pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_paired_end_read_dna_variant_calling_delly2 or \
      test_paired_end_read_dna_variant_calling_gridss or \
      test_paired_end_read_dna_variant_calling_lumpy or \
      test_paired_end_read_dna_variant_calling_manta or \
      test_paired_end_read_dna_variant_calling_mutect2_human or \
      test_paired_end_read_dna_variant_calling_mutect2_nonhuman or \
      test_paired_end_read_dna_variant_calling_sequenza or \
      test_paired_end_read_dna_variant_calling_strelka2 or \
      test_paired_end_read_dna_variant_calling_svaba"
