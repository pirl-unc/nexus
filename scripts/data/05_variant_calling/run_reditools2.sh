nexus run --nf-workflow variant_calling_reditools2.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/variant_calling_reditools2/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/05_variant_calling/samples/samples_reditools2.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr17.fa.gz \
  --reference_genes_gtf_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz \
  --params_reditools2 '""' \
  --params_reditools_annotatetable '"-s 4 -c 1,2,3 -n gencode"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/tsv/