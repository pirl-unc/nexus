docker run \
  -v "/Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/":"/input/" \
  -v "/Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/":"/output/" \
  google/deepvariant:"1.6.0" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/fasta/hg38_chr17_1-8000000.fa \
  --reads=/input/bam/hg38_tp53_tumor_long_read_dna.bam \
  --output_vcf=/output/vcf/hg38_tp53_tumor_long_read_dna_deepvariant.vcf \
  --output_gvcf=/output/vcf/hg38_tp53_tumor_long_read_dna_deepvariant.gvcf \
  --num_shards=2