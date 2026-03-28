nexus run --nf-workflow variant_calling_manta-somatic.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/variant_calling_manta/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/05_variant_calling/samples/samples_manta.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr17.fa.gz \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/ \
  --params_manta_config '""' \
  --params_manta_run '""'