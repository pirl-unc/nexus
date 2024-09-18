TRANSCRIPT_FASTA_FILE=../../../test/data/fasta/gencode_v41_tp53_transcripts.fa
OUTPUT_INDEX_FILE=../../../test/data/indices/kallisto/kallisto_gencode_v41_tp53_index
kallisto index \
  --index=$OUTPUT_INDEX_FILE \
  --threads=2 \
  $TRANSCRIPT_FASTA_FILE