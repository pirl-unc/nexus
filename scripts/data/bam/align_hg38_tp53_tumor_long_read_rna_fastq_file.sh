minimap2 \
  -ax splice:hq -uf --cs --eqx -Y -L \
  -t 2 \
  -R "@RG\\tID:hg38_tp53_tumor\\tSM:hg38_tp53_tumor\\tPL:hg38_tp53_tumor\\tLB:hg38_tp53_tumor\\tPU:hg38_tp53_tumor" \
  ../../../test/data/fasta/hg38_chr17_1-8000000.fa \
  ../../../test/data/fastq/hg38_tp53_tumor_long_read_rna.fastq.gz > ../../../test/data/bam/hg38_tp53_tumor_long_read_rna.sam

samtools sort ../../../test/data/bam/hg38_tp53_tumor_long_read_rna.sam -o ../../../test/data/bam/hg38_tp53_tumor_long_read_rna.bam
samtools index -b ../../../test/data/bam/hg38_tp53_tumor_long_read_rna.bam