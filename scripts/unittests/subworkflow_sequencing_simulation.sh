pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_sequencing_simulation_art_illumina_pe or \
      test_sequencing_simulation_nanosim_genome or \
      test_sequencing_simulation_neat or \
      test_sequencing_simulation_pbsim3_dna_1 or \
      test_sequencing_simulation_pbsim3_dna_2 or \
      test_sequencing_simulation_pbsim3_rna_1 or \
      test_sequencing_simulation_pbsim3_rna_2"
