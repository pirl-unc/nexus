nexus run --nf-workflow fastq_to_unaligned_bam.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/fastq_to_unaligned_bam/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/bam/samples_fastq_files.tsv \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/bam/ \
  --read_type single-end