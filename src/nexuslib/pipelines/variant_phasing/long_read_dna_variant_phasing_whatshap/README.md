## long_read_dna_variant_phasing_whatshap.nf

Phases variants with long-read DNA BAM files using [Whatshap](https://whatshap.readthedocs.io/en/latest/).

### Inputs / Outputs

| I/O    | Description                                                   |
|:-------|:--------------------------------------------------------------|
| Input  | `bam`, `bam.bai`, small variants `vcf` files for each sample. | 
| Output | Whatshap phased `vcf` for each sample.                        |

### Dependencies

* `Whatshap`

### Example

```
nexus run --nf-workflow long_read_dna_variant_phasing_whatshap.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --params_whatshap '"--mapq 20"'
```

### Usage

```
workflow:
    1. Run Whatshap.

usage: nexus run --nf-workflow long_read_dna_variant_phasing_whatshap.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id',
                                            'bam_file',
                                            'bam_bai_file',
                                            'small_variants_vcf_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --params_whatshap                   :   Whatshap parameters (default: '"--mapq 20"').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header                           | Description                                      |
|----------------------------------|--------------------------------------------------|
| sample_id                        | Sample ID.                                       |
| bam_file                         | Full path to `bam` file.                         |
| bam_bai_file                     | Full path to `bam.bai` file.                     |
| small_variants_vcf_file          | Full path to small variants `vcf` file.          |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--reference_genome_fasta_fai_file`
* Reference genome FASTA.FAI files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--params_whatshap`
* Refer to the [Whatshap documentation](https://whatshap.readthedocs.io/en/latest/).
* The following parameters for `Whatshap` are already included in `nexus` module for `whatshap` and should not be specified:
  * `--reference`
  * `--output`