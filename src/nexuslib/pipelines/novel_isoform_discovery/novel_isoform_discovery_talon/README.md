## novel_isoform_discovery_talon.nf

Discovers novel isoforms using [Talon](https://github.com/mortazavilab/TALON).

### Inputs / Outputs

| I/O    | Description                         |
|:-------|:------------------------------------|
| Input  | `bam` file for each sample.         | 
| Output | Talon output files for each sample. |

### Dependencies

* `Talon`

### Example

```
nexus run --nf-workflow novel_isoform_discovery_talon.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --params_talon_init_db '""' \
    --params_talon '"--create_novel_spliced_genes"'
```

### Usage

```
workflow:
    1. Run Talon.

usage: nexus run --nf-workflow novel_isoform_discovery_talon.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --reference_genes_gtf_file          :   Reference genes GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf).
    --params_talon_init_db              :   talon_initialize_database parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
    --params_talon                      :   talon parameters (default: '"--create_novel_spliced_genes"').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                  |
|--------------|------------------------------|
| sample_id    | Sample ID.                   |
| bam_file     | Full path to `bam` file.     |
| bam_bai_file | Full path to `bam.bai` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

`--reference_genome_fasta_fai_file`
* Reference genome FASTA.FAI files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

`--reference_genes_gtf_file`
* Reference GTF files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

`--params_talon_init_db`
* Refer to the [Talon documentation](https://github.com/mortazavilab/TALON).
* The following parameters for `talon` are already included in `nexus` module for `talon_initialize_database` and should not be specified:
  * `--f`
  * `--g`
  * `--a`
  * `--o`

`--params_talon`
* Refer to the [Talon documentation](https://github.com/mortazavilab/TALON).
* The following parameters for `talon` are already included in `nexus` module for `talon` and should not be specified:
  * `--f`
  * `--db`
  * `--build`
  * `--threads`
  * `--o`
