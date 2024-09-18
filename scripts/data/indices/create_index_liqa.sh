GTF_FILE=../../../test/data/gtf/gencode_v41_tp53_annotation.gtf
OUTPUT_FILE=../../../test/data/indices/liqa/liqa_gencode_v41_tp53.refgene
liqa -task refgene \
  -ref $GTF_FILE \
  -format gtf \
  -out $OUTPUT_FILE