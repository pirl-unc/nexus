pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_variant_calling_clair3rna or \
      test_variant_calling_de_souza_github or \
      test_variant_calling_longgf or \
      test_variant_calling_pbfusion"
