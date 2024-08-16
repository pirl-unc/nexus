nexus run --nf-workflow paired-end_read_rna_alignment_star.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_alignment_star/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/bam/samples_paired-end_read_rna_fastq_files_1.tsv \
  --star_index /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/star/hg38_chr6/ \
  --params_star '"--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/bam/

nexus run --nf-workflow paired-end_read_rna_alignment_star.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_alignment_star/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/bam/samples_paired-end_read_rna_fastq_files_2.tsv \
  --star_index /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/star/hg38_chr17/ \
  --params_star '"--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/bam/
