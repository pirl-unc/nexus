/home/jinseok/software/miniconda3/envs/nexus/bin/STAR \
  --runThreadN 2 \
  --runMode genomeGenerate \
  --genomeSAindexNbases 10 \
  --genomeDir ../../../test/data/indices/star/hg38_chr17/ \
  --genomeFastaFiles ../../../test/data/fasta/hg38_chr17_1-8000000.fa \
  --sjdbGTFfile ../../../test/data/gtf/gencode_v41_tp53_annotation.gtf \
  --sjdbOverhang 149
