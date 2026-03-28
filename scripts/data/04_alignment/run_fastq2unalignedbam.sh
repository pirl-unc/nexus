nexus run --nf-workflow utilities_fastq2unalignedbam.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/utilities_fastq2unalignedbam/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/04_alignment/samples/samples_fastq2unalignedbam.tsv \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/bam/ \
  --read_type single-end