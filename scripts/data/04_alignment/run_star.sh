nexus run --nf-workflow alignment_star.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/alignment_star/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/04_alignment/samples/samples_star_1.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr17.fa.gz \
  --reference_genes_gtf_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz \
  --params_star_genomegenerate '"--genomeSAindexNbases 10 --sjdbOverhang 150"' \
  --params_star '"--genomeLoad NoSharedMemory --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic --chimSegmentMin 10 --chimOutType WithinBAM SoftClip Junctions --chimMultimapNmax 20 --chimOutJunctionFormat 0"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/bam/

nexus run --nf-workflow alignment_star.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/alignment_star/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/04_alignment/samples/samples_star_2.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/GRCh38.p14.genome.chr6.fa.gz \
  --reference_genes_gtf_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr6.gtf.gz \
  --params_star '"--genomeLoad NoSharedMemory --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic --chimSegmentMin 10 --chimOutType WithinBAM SoftClip Junctions --chimMultimapNmax 20 --chimOutJunctionFormat 0"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/bam/
