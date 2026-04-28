tests=(
  test_quantification_bambu
  test_quantification_kallisto_lr
  test_quantification_kallisto_pe
  test_quantification_liqa
  test_quantification_oarfish
  test_quantification_salmon_fastq
  test_quantification_transigner
)

expr=$(printf " or %s" "${tests[@]}")
expr=${expr:4}

pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "$expr"
