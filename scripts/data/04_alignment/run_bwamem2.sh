nexus run --nf-workflow alignment_bwamem2.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/alignment_bwamem2/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/04_alignment/samples/samples_bwamem2_1.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr17.fa.gz \
  --known_sites_vcf_files '"/Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/known_sites_hg38_chr17.vcf.gz"' \
  --perform_local_indel_realignment false \
  --chromosomes '"chr17"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/bam/

nexus run --nf-workflow alignment_bwamem2.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/alignment_bwamem2/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/04_alignment/samples/samples_bwamem2_2.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr21-22.fa.gz \
  --known_sites_vcf_files '"/Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/known_sites_hg38_chr21_chr22.vcf.gz"' \
  --perform_local_indel_realignment false \
  --chromosomes '"chr21,chr22"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/bam/
