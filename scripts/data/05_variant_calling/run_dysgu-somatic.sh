nexus run --nf-workflow variant_calling_dysgu-somatic.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/variant_calling_dysgu-pe/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/05_variant_calling/samples/samples_dysgu.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr17.fa.gz \
  --params_dysgu_run '"--mode pe --min-support 3 --min-size 10 --mq 20"' \
  --params_dysgu_filter '"--support-fraction 0.05 --min-mapq 20"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/
