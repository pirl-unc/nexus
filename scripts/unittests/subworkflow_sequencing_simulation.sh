tests=(
  test_sequencing_simulation_art_illumina_pe
  test_sequencing_simulation_nanosim_genome
  test_sequencing_simulation_neat
  test_sequencing_simulation_pbsim3_dna_1
  test_sequencing_simulation_pbsim3_dna_2
  test_sequencing_simulation_pbsim3_rna_1
  test_sequencing_simulation_pbsim3_rna_2
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"
