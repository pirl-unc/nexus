## variant_calling_mutect2.nf

Identifies small somatic DNA variants (SNVs and INDELs) in paired-end DNA sequencing 
(BAM) files using [GATK4-Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2).

### Inputs / Outputs

| I/O    | Description                                                |
|:-------|:-----------------------------------------------------------|
| Input  | Tumor and normal `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.                                |

### Dependencies

* `GATK4 Mutect2`
* `Picard`

### Example

```
nexus run --nf-workflow variant_calling_mutect2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --mutect2_germline_resource_vcf_file MUTECT2_GERMLINE_RESOURCE_VCF_FILE \
    --mutect2_panel_of_normals_vcf_file MUTECT2_PANEL_OF_NORMALS_VCF_FILE \
    --getpileupsummaries_variant_vcf_file GETPILEUPSUMMARIES_VARIANT_VCF_FILE \
    --params_gatk4mutect2 '""' \
    --params_gatk4getpileupsummaries '""' \
    --chromosomes 'chr1,chr2,chr3'
```

### Usage

```
workflow:
    1. Run GATK4 Mutect2 (tumor and normal mode).
    2. Run GATK4 LearnReadOrientationModel (if --getpileupsummaries_variant_vcf_file is not empty).
    3. Run GATK4 GetPileupSummaries (if --getpileupsummaries_variant_vcf_file is not empty).
    4. Run GATK4 CalculateContamination (if --getpileupsummaries_variant_vcf_file is not empty).
    5. Run GATK4 FilterMutectCalls.
    6. Run Picard MergeVcfs,

usage: nexus run --nf-workflow variant_calling_mutect2.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file                      :   TSV file with the following columns:
                                                'sample_id', 
                                                'tumor_bam_file',
                                                'tumor_bam_bai_file',
                                                'normal_bam_file',
                                                'normal_bam_bai_file',
                                                'normal_sample_id'.
    --output_dir                            :   Directory to which output files will be symlinked.
    --reference_genome_fasta_file           :   Reference genome FASTA file.
    --mutect2_germline_resource_vcf_file    :   Germline resource VCF file. This VCF file will be supplied to gatk Mutect2 --germline-resource parameter.
    --mutect2_panel_of_normals_vcf_file     :   Panel of normals VCF file. This VCF file will be supplied to gatk Mutect2 --panel-of-normals parameter.
    --getpileupsummaries_variant_vcf_file   :   GetPileupSummaries variant VCF file.

optional arguments:
    --params_gatk4mutect2                   :   GATK4 Mutect2 parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
    --params_gatk4getpileupsummaries        :   GATK4 GetPileupSummaries parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
    --chromosomes                           :   Chromosomes to parallelize 
                                                (default: 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM').
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header              | Description                                      |
|---------------------|--------------------------------------------------|
| sample_id           | Sample ID.                                       |
| tumor_bam_file      | Full path to tumor `bam` file.                   |
| tumor_bam_bai_file  | Full path to tumor `bam.bai` file.               |
| normal_bam_file     | Full path to normal `bam` file.                  |
| normal_bam_bai_file | Full path to normal `bam.bai` file.              |
| normal_sample_id    | Normal sample ID (sample ID in normal BAM file). |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta` or `.fasta.gz`) file.

`--mutect2_germline_resource_vcf_file`
* Germline resource VCF file. 
* Recommended VCF file for `hg38` (download from [here](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/)):
  - `af-only-gnomad.hg38.vcf`

`--getpileupsummaries_variant_vcf_file`
* Panel of normals VCF file. 
* Recommended VCF file for `hg38` (download from [here](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/)):
  - `small_exac_common_3.hg38.vcf`

`--mutect2_panel_of_normals_vcf_file`
* Panel of normals VCF file. 
* Recommended VCF file for `hg38` (download from [here](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/)):
  - `1000g_pon.hg38.vcf`

`--params_gatk4mutect2`
* Refer to the [GATK4 Mutect2 documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2).
* The following parameters for `gatk` are already included in `nexus` module for `gatk mutect2` and should not be specified:
  * `-I`
  * `-normal`
  * `--native-pair-hmm-threads`
  * `--intervals`
  * `--f1r2-tar-gz`
  * `-O`

`--params_gatk4getpileupsummaries`
* Refer to the [GATK4 GetPileupSummaries documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2).
* The following parameters for `gatk` are already included in `nexus` module for `gatk mutect2` and should not be specified:
  * `-I`
  * `-V`
  * `-L`
  * `-O`
