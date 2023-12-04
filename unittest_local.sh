# Exclude test_long_read_dna_small_variants_singularity; instead test docker
# Exclude test_long_read_dna_structural_variants_github because pbsv only runs on linux
# Exclude test_paired_end_read_dna_alignment_bwa_mem2.py because abra2 only runs on linux
# Exclude test_paired_end_read_dna_somatic_small_variants_human_singularity because strelka2 only runs on linux
# Exclude test_paired_end_read_dna_somatic_structural_variants because delly2 and lumpyexpress only run on linux
pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "not test_long_read_dna_small_variants_singularity and \
      not test_long_read_dna_structural_variants_github and \
      not test_paired_end_read_dna_alignment_bwa_mem2 and \
      not test_paired_end_read_dna_somatic_small_variants_human_singularity and \
      not test_paired_end_read_dna_somatic_structural_variants"