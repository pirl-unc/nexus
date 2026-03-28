nexus run --nf-workflow variant_calling_delly2-sr-somatic.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/variant_calling_delly2/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/05_variant_calling/samples/samples_delly2.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr17.fa.gz \
  --exclude_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/delly2/human.hg38.excl.tsv \
  --params_delly2call '"--map-qual 20"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/