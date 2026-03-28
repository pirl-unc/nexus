pytest \
  -s \
  --cov-report=term-missing \
  --cov=nexuslib \
  test/ \
  -k "test_quantification_bambu or \
      test_quantification_kallisto_lr or \
      test_quantification_kallisto_pe or \
      test_quantification_liqa or \
      test_quantification_oarfish or \
      test_quantification_salmon_fastq or \
      test_quantification_transigner"
