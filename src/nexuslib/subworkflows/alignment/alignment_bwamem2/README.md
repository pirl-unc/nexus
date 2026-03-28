## alignment_bwamem2.nf

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
nexus run --nf-workflow alignment_bwamem2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --abra2_targets_bed_file ABRA2_TARGETS_BED_FILE \
    --known_sites_vcf_files VCF_FILE1,VCF_FILE2
```

### Usage

```
workflow:
    1. Index fasta and vcf files.
    2. Align paired-end reads to a reference genome using bwa-mem2.
    3. Sort sam to bam files.
    4. Perform local realignment using abra2 (optional).
    5. Add mate score tags using samtools.
    6. Mark PCR duplicates using samtools.
    7. Calculate base recalibration scores using gatk4.
    8. Apply base recalibration scores using gatk4.

usage: nexus run --nf-workflow alignment_bwamem2.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --abra2_targets_bed_file            :   ABRA2 targets BED file.
    --known_sites_vcf_files             :   GATK4 BaseRecalibrator --known-sites files. At least one VCF file must be supplied. Note that VCF files should be separated by commas.

optional arguments:
    --chromosomes                       :   Chromosomes to recalibrate using GATK4 (default:
                                            'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM').
    --platform_tag                      :   Platform tag (default: 'illumina').
    --platform_unit_tag                 :   Platform unit tag (default: 'unknown').
    --library_tag                       :   Library tag (default: 'unknown').
    --perform_local_indel_realignment   :   Perform local INDEL realignment (default: true).
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

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta`) file. An uncompressed `.fasta` file should be supplied.

`--abra2_targets_bed_file`
* A bed file with a list of `chromosome`, `start`, and `end` coordinates. 
* ABRA2 target BED files from a GENCODE GTF file can be created by running `nexus_create_abra2_targets_bed_file` (`nexus` built-in executable).

`--known_sites_vcf_files`
* A list of VCF files separated by comma.
* Recommended VCF files for hg38 (download from [here](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/)):
  - `Homo_sapiens_assembly38.dbsnp138.vcf`
  - `Homo_sapiens_assembly38.known_indels.vcf`
  - `1000G_phase1.snps.high_confidence.hg38.vcf`
  - `Mills_and_1000G_gold_standard.indels.hg38.vcf`
* Recommended VCF file for mm39 (download from [here](http://ftp.ensembl.org/pub/release-114/variation/vcf/mus_musculus/)):
  - `mus_musculus.vcf.gz`
