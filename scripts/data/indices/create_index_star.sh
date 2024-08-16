STAR \
  --runThreadN 2 \
  --runMode genomeGenerate \
  --genomeSAindexNbases 10 \
  --genomeDir ../../../test/data/indices/star/hg38_chr17/ \
  --genomeFastaFiles ../../../test/data/fasta/hg38_chr17_1-8M.fa \
  --sjdbGTFfile ../../../test/data/gtf/gencode_v41_tp53_annotation.gtf \
  --sjdbOverhang 149

STAR \
  --runThreadN 2 \
  --runMode genomeGenerate \
  --genomeSAindexNbases 10 \
  --genomeDir ../../../test/data/indices/star/hg38_chr6/ \
  --genomeFastaFiles ../../../test/data/fasta/hg38_chr6.fa \
  --sjdbGTFfile ../../../test/data/gtf/gencode_v41_hla_annotation.gtf \
  --sjdbOverhang 149
