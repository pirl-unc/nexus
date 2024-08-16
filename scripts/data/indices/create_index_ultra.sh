gff3sort.pl ../../../test/data/gtf/gencode_v41_tp53_annotation.gtf > ../../../test/data/gtf/gencode_v41_tp53_annotation_sorted.gtf
uLTRA index \
  --disable_infer \
  ../../../test/data/fasta/hg38_chr17_1-8M.fa \
  ../../../test/data/gtf/gencode_v41_tp53_annotation_sorted.gtf \
  ../../../test/data/indices/ultra/hg38_chr17/
