## long_read_dna_haplotyping_whatshap.nf

Haplotype long DNA BAM files using [Whatshap](https://whatshap.readthedocs.io/en/latest/).

### Inputs / Outputs

| I/O    | Description                                                             |
|:-------|:------------------------------------------------------------------------|
| Input  | `bam`, `bam.bai`, and `phased.vcf` files for each sample.               | 
| Output | Haplotyped `bam` and `bam.bai` files for each sample.                   |

### Dependencies

* `WhatsHap`

### Example
```
nexus run --nf-workflow long_read_dna_haplotyping_whatshap.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --params_whatshap '"--ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads 4"'
```

### Usage

```
workflow:
    1. Run Whatshap 'haplotag' command.

usage: nexus run --nf-workflow long_read_dna_haplotyping_whatshap.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id',
                                            'bam_file',
                                            'bam_bai_file',
                                            'phased_small_variants_vcf_file',
                                            'phased_small_variants_vcf_tbi_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --params_whatshap                   :   Whatshap parameters (default: '"--ignore-read-groups --tag-supplementary --skip-missing-contigs --output-threads 4"').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow `config` file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header    | Description              |
|-----------|--------------------------|
| sample_id | Sample ID.               |
| bam_file  | Full path to `bam` file. |
| bam_bai_file | Full path to `bam.bai` file. |
| phased_small_variants_vcf_file | Full path to `phased.vcf` file. |
| phased_small_variants_vcf_tbi_file | Full path to `phased.vcf.tbi` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--reference_genome_fasta_fai_file`
* Reference genome FASTA.FAI files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--params_whatshap`
* Refer to the [WhatsHap documentation](https://whatshap.readthedocs.io/en/latest/).
* The following parameters for `whatshap haplotag` are already included in `nexus` module and should not be specified:
  * `--reference`
  * `-o`
