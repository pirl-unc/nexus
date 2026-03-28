nexus run --nf-workflow variant_calling_mutect2.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/variant_calling_mutect2/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/05_variant_calling/samples/samples_mutect2.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr17.fa.gz \
  --mutect2_germline_resource_vcf_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/small_exac_common_3.hg38.vcf \
  --mutect2_panel_of_normals_vcf_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/known_sites_hg38_chr17.vcf.gz \
  --getpileupsummaries_variant_vcf_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/small_exac_common_3.hg38.vcf \
  --chromosomes '"chr17"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/