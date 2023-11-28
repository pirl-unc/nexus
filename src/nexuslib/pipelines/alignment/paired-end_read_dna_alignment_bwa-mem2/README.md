## paired_end_read_dna_alignment_bwa-mem2.nf

Aligns paired-end short DNA reads (Illumina) to a reference genome using [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2).

### Inputs / Outputs

| Input(s)                    | Output(s)                                                                                      |
| --------------------------- | ---------------------------------------------------------------------------------------------- |
| `R1 and R2 FASTQ.GZ` files  | Sorted, INDEL realigned, PCR-duplicate marked, and base recalibrated `BAM` and `BAM.BAI` files |

### Dependencies

* `bwa-mem2`
* `samtools`
* `abra2`
* `gatk4`

### Usage

```shell
nextflow run paired_end_read_dna_alignment_bwa-mem2.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir                    :   Directory to which output files will be copied.
    --reference_genome_fasta_file   :   Reference genome FASTA file.
    --bwa_mem2                      :   bwa-mem2 path.
    --samtools                      :   samtools path.
    --abra2                         :   abra2 (.jar) path.
    --abra2_targets_bed_file        :   ABRA2 targets BED file.
    --gatk4                         :   GATK4 path.
    --chromosomes                   :   Chromosomes to recalibrate using GATK4 (separated by comma; e.g. 'chr1,chr2,chr3').
    --chromosomes_count             :   Chromosomes count.
    --known_sites_vcf_files         :   Known sites (SNP and INDEL) VCF files (separated by comma).

optional arguments:
    --platform_tag                  :   Platform tag (default: 'illumina').
    --platform_unit_tag             :   Platform unit tag (default: 'unknown').
    --library_tag                   :   Library tag (default: 'unknown').
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`--sample_tsv_file`

| Header       | Description                     |
| ------------ | ------------------------------- |
| sample_id    | Sample ID                       |
| fastq_file_1 | Full path to R1 `FASTQ.GZ` file |
| fastq_file_2 | Full path to R2 `FASTQ.GZ` file |

`--abra2_targets_bed_file`

A bed file with a list of `chromosome`, `start`, and `end` coordinates.

`--known_sites_vcf_files`

Supply the following VCF files for `hg38`:
- [Homo_sapiens_assembly38.known_indels.vcf](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false)
- [Mills_and_1000G_gold_standard.indels.hg38.vcf](ftp://ftp.broadinstitute.org/bundle/hg38/)
- [dbsnp_*.hg38.vcf](ftp://ftp.broadinstitute.org/bundle/hg38/)

