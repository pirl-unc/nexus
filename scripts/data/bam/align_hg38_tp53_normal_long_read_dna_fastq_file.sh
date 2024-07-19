minimap2 \
  -ax map-hifi --cs --eqx -Y -L \
  -t 2 \
  -R "@RG\\tID:hg38_tp53_normal\\tSM:hg38_tp53_normal\\tPL:hg38_tp53_normal\\tLB:hg38_tp53_normal\\tPU:hg38_tp53_normal" \
  ../../../test/data/fasta/hg38_chr17_1-8000000.fa \
  ../../../test/data/fastq/hg38_tp53_normal_long_read_dna.fastq.gz > ../../../test/data/bam/hg38_tp53_normal_long_read_dna.sam

samtools sort ../../../test/data/bam/hg38_tp53_normal_long_read_dna.sam -o ../../../test/data/bam/hg38_tp53_normal_long_read_dna.bam
samtools index -b ../../../test/data/bam/hg38_tp53_normal_long_read_dna.bam