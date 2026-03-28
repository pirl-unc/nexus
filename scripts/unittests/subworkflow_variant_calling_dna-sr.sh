pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_variant_calling_clair3_2 or \
      test_variant_calling_clairs_2 or \
      test_variant_calling_deepsomatic_github_2 or \
      test_variant_calling_deepvariant_github_2 or \
      test_variant_calling_delly2_sr_germline or \
      test_variant_calling_delly2_sr_somatic or \
      test_variant_calling_dysgu_germline_pe or \
      test_variant_calling_dysgu_somatic_pe or \
      test_variant_calling_gridss2_germline or \
      test_variant_calling_gridss2_somatic or \
      test_variant_calling_haplotypecaller or \
      test_variant_calling_lumpy_germline or \
      test_variant_calling_lumpy_somatic or \
      test_variant_calling_manta_germline or \
      test_variant_calling_manta_somatic or \
      test_variant_calling_mutect2_1 or \
      test_variant_calling_mutect2_2 or \
      test_variant_calling_octopus_germline or \
      test_variant_calling_octopus_somatic or \
      test_variant_calling_pindel or \
      test_variant_calling_sequenza or \
      test_variant_calling_strelka2_germline or \
      test_variant_calling_strelka2_somatic or \
      test_variant_calling_svaba"
