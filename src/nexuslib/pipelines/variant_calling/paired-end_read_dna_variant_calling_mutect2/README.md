## paired_end_read_dna_variant_calling_mutect2.nf

Identifies small somatic DNA variants (SNVs and INDELs) in paired-end DNA sequencing 
(BAM) files using [GATK4-Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2).

### Inputs / Outputs

| I/O    | Description                                   |
|:-------|:----------------------------------------------|
| Input  | Tumor and normal `bam` files for each sample. | 
| Output | `vcf` file for each sample.                   |

### Dependencies

* `GATK4 Mutect2`
* `Picard`

### Example

```
nexus run --nf-workflow paired-end_read_dna_variant_calling_mutect2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --reference_genome_fasta_dict_file REFERENCE_GENOME_FASTA_DICT_FILE \
    --is_human true \
    --params_gatk4mutect2 '"--germline-resource VCF_FILE --panel-of-normals VCF_FILE"' \
    --params_gatk4getpileupsummaries '"-V VCF_FILE -L VCF_FILE"' \
    --chromosomes 'chr1,chr2,chr3'
```

### Usage

```
workflow:
    human:
        1.  Run GATK4 Mutect2 (tumor and normal mode).
            Run GATK4 LearnReadOrientationModel.
            Run GATK4 GetPileupSummaries.
            Run GATK4 CalculateContamination.
            Run GATK4 FilterMutectCalls.
            Run Picard MergeVcfs,
    non-human:
        1.  Run GATK4 Mutect2 (tumor and normal mode).
            Run GATK4 FilterMutectCalls.
            Run Picard MergeVcfs,

usage: nexus run --nf-workflow paired-end_read_dna_variant_calling_mutect2.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file', 'tumor_sample_id', normal_sample_id'
    --output_dir                        :   Directory to which output files will be symlinked.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --reference_genome_fasta_dict_file  :   Reference genome DICT file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.dict).
    --is_human                          :   true if the samples are human. false otherwise (default: true).
    --params_gatk4mutect2               :   GATK4 Mutect2 parameters (default:
                                            '"--germline-resource /datastore/lbcfs/collaborations/pirl/seqdata/references/af-only-gnomad.hg38.vcf
                                              --panel-of-normals /datastore/lbcfs/collaborations/pirl/seqdata/references/1000g_pon.hg38.vcf"').
                                            Note that the parameters need to be wrapped in quotes.
    --params_gatk4getpileupsummaries    :   GATK4 GetPileupSummaries parameters (default:
                                            '"-V /datastore/lbcfs/collaborations/pirl/seqdata/references/small_exac_common_3.hg38.vcf
                                              -L /datastore/lbcfs/collaborations/pirl/seqdata/references/small_exac_common_3.hg38.vcf"').
                                            Note that the parameters need to be wrapped in quotes.
    --chromosomes                       :   Chromosomes to parallelize (default: 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM').
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`

| Header              | Description                           |
|---------------------|---------------------------------------|
| sample_id           | Sample ID.                            |
| tumor_bam_file      | Full path to tumor `bam` file.        |
| tumor_bam_bai_file  | Full path to tumor `bam.bai` file.    |
| normal_bam_file     | Full path to normal `bam` file.       |
| normal_bam_bai_file | Full path to normal `bam.bai` file.   |
| normal_sample_id    | Normal sample ID.                     |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--reference_genome_fasta_fai_file`
* Reference genome FASTA.FAI files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--reference_genome_fasta_dict_file`
* Reference genome FASTA.DICT files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--params_gatk4mutect2`
* Recommended values for `hg38`:
  * `--germline-resource af-only-gnomad.hg38.vcf --panel-of-normals 1000g_pon.hg38.vcf`

`--params_gatk4getpileupsummaries`
* Recommended values for `hg38`:
  * `-V small_exac_common_3.hg38.vcf -L small_exac_common_3.hg38.vcf`