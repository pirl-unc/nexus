pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_variant_calling_clair3_1 or \
      test_variant_calling_clairs_1 or \
      test_variant_calling_cutesv or \
      test_variant_calling_deepsomatic_github_1 or \
      test_variant_calling_deepvariant_github_1 or \
      test_variant_calling_delly2_lr_germline or \
      test_variant_calling_delly2_lr_somatic or \
      test_variant_calling_dysgu_somatic_lr or \
      test_variant_calling_hificnv or \
      test_variant_calling_longshot or \
      test_variant_calling_nanocaller_1 or \
      test_variant_calling_nanocaller_2 or \
      test_variant_calling_nanomonsv or \
      test_variant_calling_nanovar_1 or \
      test_variant_calling_nanovar_2 or \
      test_variant_calling_pbsv or \
      test_variant_calling_savana or \
      test_variant_calling_severus or \
      test_variant_calling_sniffles2 or \
      test_variant_calling_svisionpro"
