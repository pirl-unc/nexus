# Exclude test_paired_end_read_dna_alignment_bwa_mem2.py because abra2 only runs on linux
# Exclude test_long_read_dna_small_variants_singularity on local machine (test docker)
pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "not test_paired_end_read_dna_alignment_bwa_mem2 and not test_long_read_dna_small_variants_singularity"