/home/jinseok/software/miniconda3/envs/nexus/bin/STAR \
  --runThreadN 2 \
  --readFilesIn ../../../test/data/fastq/hg38_hla-abc_paired_end_read_rna.r1.fastq.gz ../../../test/data/fastq/hg38_hla-abc_paired_end_read_rna.r2.fastq.gz \
  --genomeDir /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/star/hg38_149bp_overhang/ \
  --outFileNamePrefix ../../../test/data/bam/hg38_hla-abc_star_ \
  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic
samtools index -b ../../../test/data/bam/hg38_hla-abc_star_Aligned.sortedByCoord.out.bam