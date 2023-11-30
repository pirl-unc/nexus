## paired_end_read_dna_alignment_bwa-mem2.nf

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

### Usage

```
workflow:
    1. Align paired-end reads to a reference genome using bwa-mem2.
    2. Sort sam to bam files.
    3. Perform local realignment using abra2.
    4. Add mate score tags using samtools.
    5. Mark PCR duplicates using samtools.
    6. Calculate base recalibration scores using gatk4.
    7. Apply base recalibration scores using gatk4.

usage: nexus run --nf-workflow paired_end_read_dna_alignment_bwa-mem2.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir                    :   Directory to which output files will be copied.
    --reference_genome_fasta_file   :   Reference genome FASTA file.
    --bwa_mem2                      :   bwa-mem2 path.
    --samtools                      :   samtools path.
    --abra2                         :   abra2 (.jar) path.
    --abra2_targets_bed_file        :   abra2 targets BED file.
    --gatk4                         :   gatk4 path.
    --gatk4_baserecalibrator_params :   gatk4 BaseRecalibrator parameters (e.g. '"--known-sites /path/dbsnp.vcf --known-sites /path/common.vcf "').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
    --chromosomes                   :   Chromosomes to recalibrate using GATK4 (separated by comma; e.g. 'chr1,chr2,chr3').

optional arguments:
    --platform_tag                  :   Platform tag (default: 'illumina').
    --platform_unit_tag             :   Platform unit tag (default: 'unknown').
    --library_tag                   :   Library tag (default: 'unknown').
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/configs)

`--sample_tsv_file`

| Header       | Description                     |
| ------------ |---------------------------------|
| sample_id    | Sample ID.                      |
| fastq_file_1 | Full path to R1 `fastq.gz` file |
| fastq_file_2 | Full path to R2 `fastq.gz` file |

`--abra2_targets_bed_file`
* A bed file with a list of `chromosome`, `start`, and `end` coordinates. 
* BED files are available in `/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/abra2/` on the LBG.

`--gatk4_baserecalibrator_params`
* VCF files for `--known-sites` are available in `/datastore/lbcfs/collaborations/pirl/seqdata/references` on the LBG.
* Recommended VCF files for hg38:
  - `dbsnp_146.hg38.vcf`
  - `1000G_phase1.snps.high_confidence.hg38.vcf`
  - `Mills_and_1000G_gold_standard.indels.hg38.vcf`
  - `Homo_sapiens_assembly38.known_indels.vcf`