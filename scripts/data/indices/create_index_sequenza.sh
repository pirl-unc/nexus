FASTA_FILE=../../../test/data/fasta/hg38_chr21_chr22.fa.gz
WINDOW_SIZE=50
OUTPUT_FILE=../../../test/data/indices/sequenza/hg38_chr21_chr22.gc50.wig.gz
sequenza-utils gc_wiggle \
  -w $WINDOW_SIZE \
  -f $FASTA_FILE \
  -o $OUTPUT_FILE
