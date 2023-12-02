# Exclude test_paired_end_read_dna_alignment_bwa_mem2.py because abra2 only runs on linux
# Exclude test_long_read_dna_small_variants_singularity; instead test docker)
# Exclude test_paired_end_read_dna_somatic_small_variants_human_docker because strelka2 only runs on linux
# Exclude test_paired_end_read_dna_somatic_small_variants_human_singularity because strelka2 only runs on linux
pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "not test_paired_end_read_dna_alignment_bwa_mem2 and \
      not test_long_read_dna_small_variants_singularity and \
      not test_paired_end_read_dna_somatic_small_variants_human_docker and \
      not test_paired_end_read_dna_somatic_small_variants_human_singularity"