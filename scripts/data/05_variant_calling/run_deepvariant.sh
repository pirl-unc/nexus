nexus run --nf-workflow variant_calling_deepvariant.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/variant_calling_deepvariant/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/05_variant_calling/samples/samples_deepvariant_1.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr17.fa.gz \
  --deepvariant_containerization docker \
  --deepvariant_input_path /Users/leework/Documents/Research/projects/ \
  --deepvariant_output_path /var/ \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

nexus run --nf-workflow variant_calling_deepvariant.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_deepvariant/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/05_variant_calling/samples/samples_deepvariant_2.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr21-22.fa.gz \
  --deepvariant_containerization docker \
  --deepvariant_input_path /Users/leework/Documents/Research/projects/ \
  --deepvariant_output_path /var/ \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/