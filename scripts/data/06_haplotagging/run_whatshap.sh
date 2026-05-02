nexus run --nf-workflow haplotagging_whatshap.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/haplotagging_whatshap/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/06_haplotagging/samples/samples_whatshap.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr17.fa.gz \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/
