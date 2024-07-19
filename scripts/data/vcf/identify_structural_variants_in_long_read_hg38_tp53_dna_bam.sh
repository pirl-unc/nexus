nexus run --nf-workflow long_read_dna_variant_calling_pbsv.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_pbsv/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/hg38_tp53_long_read_dna_samples.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8000000.fa \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

nexus run --nf-workflow long_read_dna_variant_calling_savana.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_savana/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/hg38_tp53_long_read_dna_tumor_normal_bam_files.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8000000.fa \
  --params_savana '"--length 30 --mapq 20 --depth 3"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

nexus run --nf-workflow long_read_dna_variant_calling_svisionpro.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_svisionpro/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/hg38_tp53_long_read_dna_tumor_normal_bam_files.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8000000.fa \
  --params_svisionpro '"--detect_mode somatic --preset hifi --min_supp 3 --min_mapq 20 --min_sv_size 30 --max_sv_size 1000000 --device cpu"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

nexus run --nf-workflow long_read_dna_variant_calling_severus.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_severus/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/hg38_tp53_long_read_dna_tumor_normal_bam_phased_vcf_files.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8000000.fa \
  --vntr_bed_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/severus/human_GRCh38_no_alt_analysis_set.trf.bed \
  --params_severus '"--min-support 3 --min-sv-size 30 --min-mapq 20 --output-read-ids --bp-cluster-size 50"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

