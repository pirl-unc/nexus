nexus run --nf-workflow long_read_dna_variant_phasing_whatshap.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_phasing_deepvariant/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/hg38_tp53_long_read_dna_bam_vcf_files.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8000000.fa \
  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8000000.fa.fai \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/