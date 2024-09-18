TRANSCRIPT_FASTA_FILE=$(realpath ../../../test/data/fasta/gencode_v41_tp53_transcripts.fa)
OUTPUT_INDEX_DIR=$(realpath ../../../test/data/indices/salmon/gencode_v41_tp53/)
salmon index \
  --transcripts $TRANSCRIPT_FASTA_FILE \
  --index $OUTPUT_INDEX_DIR \
  --threads 2 \
  --gencode
