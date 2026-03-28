######## DNA ########
nexus run --nf-workflow alignment_minimap2.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/alignment_minimap2/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/04_alignment/samples/samples_minimap2_1.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr17.fa.gz \
  --params_minimap2 '"-ax map-hifi --cs --eqx -Y -L"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/bam/

nexus run --nf-workflow alignment_minimap2.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/alignment_minimap2/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/04_alignment/samples/samples_minimap2_2.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr21-22.fa.gz \
  --params_minimap2 '"-ax map-hifi --cs --eqx -Y -L"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/bam/


######## RNA ########
nexus run --nf-workflow alignment_minimap2.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/alignment_minimap2/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/04_alignment/samples/samples_minimap2_3.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr17.fa.gz \
  --params_minimap2 '"-ax splice:hq -uf --cs --eqx -Y -L"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/bam/

nexus run --nf-workflow alignment_minimap2.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/alignment_minimap2/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/04_alignment/samples/samples_minimap2_4.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr6.fa.gz \
  --params_minimap2 '"-ax splice:hq -uf --cs --eqx -Y -L"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/bam/
