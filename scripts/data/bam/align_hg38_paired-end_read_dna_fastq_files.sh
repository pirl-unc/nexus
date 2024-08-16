nexus run --nf-workflow paired-end_read_dna_alignment_bwa-mem2.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_alignment_bwa-mem2/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/bam/samples_paired-end_read_dna_fastq_files_1.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa.gz \
  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa.gz.fai \
  --reference_genome_fasta_gzi_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa.gz.gzi \
  --reference_genome_fasta_dict_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.dict \
  --reference_genome_fasta_0123_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa.gz.0123 \
  --reference_genome_fasta_amb_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa.gz.amb \
  --reference_genome_fasta_ann_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa.gz.ann \
  --reference_genome_fasta_bwt_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa.gz.bwt.2bit.64 \
  --reference_genome_fasta_pac_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa.gz.pac \
  --known_sites_files /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/known_sites_hg38_chr17.vcf \
  --perform_local_indel_realignment false \
  --chromosomes '"chr17"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/bam/

#bwa-mem2 mem -t 4 \
#    -R "@RG\tID:sample105tumor\tSM:sample105tumor\tPL:ILMN\tLB:unknown\tPU:ILMN" \
#    ../../../test/data/fasta/hg38_chr21_chr22.fa.gz \
#    ../../../test/data/fastq/sample105tumor_paired-end_read_dna.r1.fastq.gz \
#    ../../../test/data/fastq/sample105tumor_paired-end_read_dna.r2.fastq.gz > ../../../test/data/bam/sample105tumor.sam
#samtools sort -@ 4 -m 1G -O bam -o ../../../test/data/bam/sample105tumor.bam ../../../test/data/bam/sample105tumor.sam
#samtools index -b ../../../test/data/bam/sample105tumor.bam ../../../test/data/bam/sample105tumor.bam.bai
#rm ../../../test/data/bam/sample105tumor.sam
#
#bwa-mem2 mem -t 4 \
#    -R "@RG\tID:sample105normal\tSM:sample105normal\tPL:ILMN\tLB:unknown\tPU:ILMN" \
#    ../../../test/data/fasta/hg38_chr21_chr22.fa.gz \
#    ../../../test/data/fastq/sample105normal_paired-end_read_dna.r1.fastq.gz \
#    ../../../test/data/fastq/sample105normal_paired-end_read_dna.r2.fastq.gz > ../../../test/data/bam/sample105normal.sam
#samtools sort -@ 4 -m 1G -O bam -o ../../../test/data/bam/sample105normal.bam ../../../test/data/bam/sample105normal.sam
#samtools index -b ../../../test/data/bam/sample105normal.bam ../../../test/data/bam/sample105normal.bam.bai
#rm ../../../test/data/bam/sample105normal.sam