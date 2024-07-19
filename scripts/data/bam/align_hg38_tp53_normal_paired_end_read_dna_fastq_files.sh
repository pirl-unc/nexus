FASTQ_FILE_1=../../../test/data/fastq/hg38_tp53_normal_paired-end_read_dna.r1.fastq.gz
FASTQ_FILE_2=../../../test/data/fastq/hg38_tp53_normal_paired-end_read_dna.r2.fastq.gz
FASTA_FILE=../../../test/data/fasta/hg38_chr17_1-8000000.fa
SAM_FILE=../../../test/data/bam/hg38_tp53_normal_paired-end_read_dna.sam
BAM_FILE=../../../test/data/bam/hg38_tp53_normal_paired-end_read_dna.bam
BAM_BAI_FILE=../../../test/data/bam/hg38_tp53_normal_paired-end_read_dna.bam.bai

samtools faidx $FASTA_FILE
bwa-mem2 index $FASTA_FILE
gatk CreateSequenceDictionary -R $FASTA_FILE

bwa-mem2 mem -t 4 \
    -R "@RG\tID:synthetic001\tSM:sample001normal\tPL:ILMN\tLB:unknown\tPU:ILMN" \
    $FASTA_FILE \
    $FASTQ_FILE_1 $FASTQ_FILE_2 > $SAM_FILE

samtools sort -@ 4 -m 1G -O bam -o $BAM_FILE $SAM_FILE
samtools index -b $BAM_FILE $BAM_BAI_FILE

rm $SAM_FILE