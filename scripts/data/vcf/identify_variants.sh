# Long-read Variant Calling
## ClairS
#nexus run --nf-workflow long_read_dna_variant_calling_clairs.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_clairs/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_long_read_dna_tumor_normal_bam_files.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M.fa.gz \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M.fa.gz.fai \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

## CuteSV
#nexus run --nf-workflow long_read_dna_variant_calling_cutesv.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_cutesv/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_long_read_dna_bam_files.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.gz \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.gz.fai \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

# DeepVariant
nexus run --nf-workflow long_read_dna_variant_calling_deepvariant.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_deepvariant/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_long_read_dna_bam_files.tsv \
  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa \
  --containerization docker \
  --deepvariant_input_path /Users/leework/Documents/Research/projects/ \
  --deepvariant_output_path /var/ \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

## Pbsv
#nexus run --nf-workflow long_read_dna_variant_calling_pbsv.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_pbsv/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_long_read_dna_bam_files.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

## Savana
#nexus run --nf-workflow long_read_dna_variant_calling_savana.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_savana/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_long_read_dna_tumor_normal_bam_files.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.gz \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.gz.fai \
#  --custom_params_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/savana/savana_classification_parameters.json \
#  --params_savana_run '"--length 30 --mapq 20 --min_support 3"' \
#  --params_savana_classify '""' \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

## Severus
#nexus run --nf-workflow long_read_dna_variant_calling_severus.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_severus/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_long_read_dna_tumor_normal_bam_phased_vcf_files.tsv \
#  --vntr_bed_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/severus/human_GRCh38_no_alt_analysis_set.trf.bed \
#  --params_severus '"--min-support 3 --min-sv-size 30 --min-mapq 20 --output-read-ids --bp-cluster-size 50"' \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

## Sniffles2
#nexus run --nf-workflow long_read_dna_variant_calling_sniffles2.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_sniffles2/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_long_read_dna_bam_files.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.gz \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.gz.fai \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

## SVIM
#nexus run --nf-workflow long_read_dna_variant_calling_svim.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_svim/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_long_read_dna_bam_files.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.gz \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.gz.fai \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

## SVision-Pro
#nexus run --nf-workflow long_read_dna_variant_calling_svisionpro.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/long_read_dna_variant_calling_svisionpro/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_long_read_dna_tumor_normal_bam_files.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.gz \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.gz.fai \
#  --svisionpro_model_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/svisionpro/model_liteunet_256_8_16_32_32_32.pth \
#  --params_svisionpro '"--detect_mode somatic --preset hifi --min_supp 3 --min_mapq 20 --min_sv_size 30 --max_sv_size 1000000 --device cpu"' \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

# Paired-end read Variant Calling
#nexus run --nf-workflow paired-end_read_dna_variant_calling_delly2.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_variant_calling_delly2/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_paired-end_read_dna_tumor_normal_bam_files_1.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.gz \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M_chr18_1-9M_hpv16.fa.gz.fai \
#  --exclude_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/delly2/human.hg38.excl.tsv \
#  --params_delly2call '"--map-qual 20"' \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/
#nexus run --nf-workflow paired-end_read_dna_variant_calling_delly2.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_variant_calling_delly2/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_paired-end_read_dna_tumor_normal_bam_files_2.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr21_chr22.fa.gz \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr21_chr22.fa.gz.fai \
#  --exclude_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/delly2/human.hg38.excl.tsv \
#  --params_delly2call '"--map-qual 20"' \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

#nexus run --nf-workflow paired-end_read_dna_variant_calling_gridss.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_variant_calling_gridss/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_paired-end_read_dna_tumor_normal_bam_files_1.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa.fai \
#  --params_gridss '""' \
#  --params_gridss_somatic_filter '"--ref BSgenome.Hsapiens.UCSC.hg38"' \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/
#nexus run --nf-workflow paired-end_read_dna_variant_calling_gridss.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_variant_calling_gridss/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_paired-end_read_dna_tumor_normal_bam_files_2.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr21_chr22.fa \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr21_chr22.fa.fai \
#  --params_gridss '""' \
#  --params_gridss_somatic_filter '"--ref BSgenome.Hsapiens.UCSC.hg38"' \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

#nexus run --nf-workflow paired-end_read_dna_variant_calling_lumpy.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_variant_calling_lumpy/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_paired-end_read_dna_tumor_normal_bam_files_1.tsv \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/
#nexus run --nf-workflow paired-end_read_dna_variant_calling_lumpy.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_variant_calling_lumpy/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_paired-end_read_dna_tumor_normal_bam_files_2.tsv \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

#nexus run --nf-workflow paired-end_read_dna_variant_calling_manta.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_variant_calling_manta/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_paired-end_read_dna_tumor_normal_bam_files_1.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa.fai \
#  --params_manta_config '""' \
#  --params_manta_run '""' \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/
#nexus run --nf-workflow paired-end_read_dna_variant_calling_manta.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_variant_calling_manta/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_paired-end_read_dna_tumor_normal_bam_files_2.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr21_chr22.fa \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr21_chr22.fa.fai \
#  --params_manta_config '""' \
#  --params_manta_run '""' \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

#nexus run --nf-workflow paired-end_read_dna_variant_calling_mutect2.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_variant_calling_mutect2/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_paired-end_read_dna_tumor_normal_bam_files_1.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa.fai \
#  --reference_genome_fasta_dict_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.dict \
#  --is_human 0 \
#  --params_gatk4mutect2 '""' \
#  --params_gatk4getpileupsummaries '""' \
#  --chromosomes 'chr17' \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/
#nexus run --nf-workflow paired-end_read_dna_variant_calling_mutect2.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_variant_calling_mutect2/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_paired-end_read_dna_tumor_normal_bam_files_2.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr21_chr22.fa \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr21_chr22.fa.fai \
#  --reference_genome_fasta_dict_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr21_chr22.dict \
#  --is_human 0 \
#  --params_gatk4mutect2 '""' \
#  --params_gatk4getpileupsummaries '""' \
#  --chromosomes '"chr21,chr22"' \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

#nexus run --nf-workflow paired-end_read_dna_variant_calling_sequenza.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_variant_calling_sequenza/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_paired-end_read_dna_tumor_normal_bam_files_2.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr21_chr22.fa \
#  --reference_genome_wig_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/sequenza/hg38_chr21_chr22.gc50.wig.gz \
#  --assembly hg38 \
#  --chromosomes '"chr21 chr22"' \
#  --params_sequenza_bam2seqz '"-N 10 --qformat sanger"' \
#  --params_sequenza_seqzbinning '"--window 50"' \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

#nexus run --nf-workflow paired-end_read_dna_variant_calling_strelka2.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_variant_calling_strelka2/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_paired-end_read_dna_tumor_normal_bam_files_1.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/fasta/hg38_chr17_1-8M.fa.fai \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/

#nexus run --nf-workflow paired-end_read_dna_variant_calling_svaba.nf \
#  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
#  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/paired-end_read_dna_variant_calling_svaba/ \
#  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/vcf/samples_hg38_paired-end_read_dna_tumor_normal_bam_files_1.tsv \
#  --reference_genome_fasta_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/svaba/hg38_chr17_1-8M.fa \
#  --reference_genome_fasta_amb_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/svaba/hg38_chr17_1-8M.fa.amb \
#  --reference_genome_fasta_ann_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/svaba/hg38_chr17_1-8M.fa.ann \
#  --reference_genome_fasta_bwt_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/svaba/hg38_chr17_1-8M.fa.bwt \
#  --reference_genome_fasta_fai_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/svaba/hg38_chr17_1-8M.fa.fai \
#  --reference_genome_fasta_pac_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/svaba/hg38_chr17_1-8M.fa.pac \
#  --reference_genome_fasta_sa_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/indices/svaba/hg38_chr17_1-8M.fa.sa \
#  --params_svaba '"--hp --read-tracking"' \
#  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/vcf/
