nexus run --nf-workflow variant_calling_circexplorer2.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/variant_calling_circexplorer2/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/05_variant_calling/samples/samples_circexplorer2.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr17.fa.gz \
  --circexplorer2_gene_annotation_txt_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/circexplorer2/hg38_kg/hg38_kg.txt \
  --params_circexplorer2_parse '"-t STAR"' \
  --params_circexplorer2_annotate '""' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/tsv/
