## paired-end_read_dna_alignment_bwa-mem2.nf

Aligns paired-end short DNA reads (Illumina) to a reference genome using [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2).

### Inputs / Outputs

| I/O    | Description                                                                                                                     |
|:-------|:--------------------------------------------------------------------------------------------------------------------------------|
| Input  | `R1 and R2 fastq.gz` files for each sample.                                                                                     | 
| Output | Sorted, INDEL-realigned, PCR_duplicate marked, and recalibrated `bam` and `bam.bai` files for each sample. |

### Dependencies

* `bwa-mem2`
* `samtools`
* `abra2`
* `gatk4`

### Example

```
nexus run --nf-workflow paired-end_read_dna_alignment_bwa-mem2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --reference_genome_fasta_gzi_file REFERENCE_GENOME_FASTA_GZI_FILE \
    --reference_genome_fasta_dict_file REFERENCE_GENOME_FASTA_DICT_FILE \
    --reference_genome_fasta_0123_file REFERENCE_GENOME_FASTA_0123_FILE \
    --reference_genome_fasta_amb_file REFERENCE_GENOME_FASTA_AMB_FILE \
    --reference_genome_fasta_ann_file REFERENCE_GENOME_FASTA_ANN_FILE \
    --reference_genome_fasta_bwt_file REFERENCE_GENOME_FASTA_BWT_FILE \
    --reference_genome_fasta_pac_file REFERENCE_GENOME_FASTA_PAC_FILE \
    --abra2_targets_bed_file ABRA2_TARGETS_BED_FILE \
    --known_sites_files VCF_FILE1,VCF_FILE2 \
    --chromosomes 'chr1,chr2,chr3'
```

### Usage

```
workflow:
    1. Align paired-end reads to a reference genome using bwa-mem2.
    2. Sort sam to bam files.
    3. (Optional) Perform local realignment using abra2.
    4. Add mate score tags using samtools.
    5. Mark PCR duplicates using samtools.
    6. Calculate base recalibration scores using gatk4.
    7. Apply base recalibration scores using gatk4.

usage: nexus run --nf-workflow paired-end_read_dna_alignment_bwa-mem2.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --reference_genome_fasta_gzi_file   :   Reference genome FASTA.GZ.GZI file (default: '').
    --reference_genome_fasta_dict_file  :   Reference genome DICT file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.dict).
    --reference_genome_fasta_0123_file  :   Reference genome FASTA.0123 file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.0123).
    --reference_genome_fasta_amb_file   :   Reference genome FASTA.AMB file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.amb).
    --reference_genome_fasta_ann_file   :   Reference genome FASTA.ANN file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.ann).
    --reference_genome_fasta_bwt_file   :   Reference genome FASTA.BWT file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.bwt.2bit.64).
    --reference_genome_fasta_pac_file   :   Reference genome FASTA.PAC file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.pac).
    --abra2_targets_bed_file            :   ABRA2 targets BED file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/abra2/gencode-v41-annotation-abra2-exon-targets.bed).
    --known_sites_files                 :   GATK4 BaseRecalibrator --known-sites files
                                            (default: '/datastore/lbcfs/collaborations/pirl/seqdata/references/dbsnp_146.hg38.vcf,/datastore/lbcfs/collaborations/pirl/seqdata/references/1000G_phase1.snps.high_confidence.hg38.vcf,/datastore/lbcfs/collaborations/pirl/seqdata/references/Mills_and_1000G_gold_standard.indels.hg38.vcf,/datastore/lbcfs/collaborations/pirl/seqdata/references/Homo_sapiens_assembly38.known_indels.vcf').
                                            Note that files should be separated by commas.
    --chromosomes                       :   Chromosomes to recalibrate using GATK4 (default:
                                            'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM').
    --platform_tag                      :   Platform tag (default: 'illumina').
    --platform_unit_tag                 :   Platform unit tag (default: 'unknown').
    --library_tag                       :   Library tag (default: 'unknown').
    --perform_local_indel_realignment   :   Perform local INDEL realignment (default: true).
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                     |
| ------------ |---------------------------------|
| sample_id    | Sample ID.                      |
| fastq_file_1 | Full path to R1 `fastq.gz` file |
| fastq_file_2 | Full path to R2 `fastq.gz` file |

`--abra2_targets_bed_file`
* A bed file with a list of `chromosome`, `start`, and `end` coordinates. 
* BED files are available in `/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/abra2/` on LBG.

`--known_sites_files`
* VCF files are available in `/datastore/lbcfs/collaborations/pirl/seqdata/references` on the LBG.
* Recommended VCF files for hg38:
  - `dbsnp_146.hg38.vcf`
  - `1000G_phase1.snps.high_confidence.hg38.vcf`
  - `Mills_and_1000G_gold_standard.indels.hg38.vcf`
  - `Homo_sapiens_assembly38.known_indels.vcf`