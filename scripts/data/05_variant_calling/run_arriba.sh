nexus run --nf-workflow variant_calling_arriba.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/variant_calling_arriba/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/05_variant_calling/samples/samples_arriba.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr17.fa \
  --reference_genes_gtf_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz \
  --protein_domains_gff3_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3 \
  --params_arriba '"-S 3 -f blacklist -i chr*"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/tsv/