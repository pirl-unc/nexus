FASTA_FILE=../../../test/data/fasta/hg38_chr21-22.fa.gz
WINDOW_SIZE=50
OUTPUT_FILE=../../../test/data/indices/sequenza/hg38_chr21-22.gc50.wig.gz
sequenza-utils gc_wiggle \
  -w $WINDOW_SIZE \
  -f $FASTA_FILE \
  -o $OUTPUT_FILE

FASTA_FILE=../../../test/data/fasta/hg38_chr17_1-8000000.fa
WINDOW_SIZE=10000
OUTPUT_FILE=../../../test/data/indices/sequenza/hg38_chr17_1-8000000.gc10000.wig.gz
sequenza-utils gc_wiggle \
  -w $WINDOW_SIZE \
  -f $FASTA_FILE \
  -o $OUTPUT_FILE
