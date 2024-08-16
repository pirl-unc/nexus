FASTA_FILE=../../../test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.gz
samtools faidx $FASTA_FILE
bwa-mem2 index $FASTA_FILE
gatk CreateSequenceDictionary -R $FASTA_FILE

FASTA_FILE=../../../test/data/fasta/hg38_chr17_1-8M_chr18_1-9M.fa.gz
samtools faidx $FASTA_FILE
bwa-mem2 index $FASTA_FILE
gatk CreateSequenceDictionary -R $FASTA_FILE

FASTA_FILE=../../../test/data/fasta/hg38_chr17_1-8M.fa.gz
samtools faidx $FASTA_FILE
bwa-mem2 index $FASTA_FILE
gatk CreateSequenceDictionary -R $FASTA_FILE

FASTA_FILE=../../../test/data/fasta/hg38_chr17_1-8M.fa
samtools faidx $FASTA_FILE

FASTA_FILE=../../../test/data/fasta/hg38_chr21_chr22.fa.gz
samtools faidx $FASTA_FILE
bwa-mem2 index $FASTA_FILE
gatk CreateSequenceDictionary -R $FASTA_FILE

FASTA_FILE=../../../test/data/fasta/hg38_chr21_chr22.fa
samtools faidx $FASTA_FILE

