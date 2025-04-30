## novel_isoform_discovery_sqanti3_fastq_mode.nf

Discovers novel isoforms using [Sqanti3](https://github.com/ConesaLab/SQANTI3).

### Inputs / Outputs

| I/O    | Description                           |
|:-------|:--------------------------------------|
| Input  | `fastq.gz` file for each sample.      | 
| Output | Sqanti3 output files for each sample. |

### Dependencies

* `Sqanti3`

### Example

```
nexus run --nf-workflow novel_isoform_discovery_sqanti3_fastq_mode.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --params_sqanti3_qc '"--aligner_choice minimap2 --report html --force_id_ignore --isoAnnotLite"' \
    --params_sqanti3_filter '"ml"'
```

### Usage

```
workflow:
    1. Run SQANTI3 (FASTQ mode).

usage: nexus run --nf-workflow novel_isoform_discovery_sqanti3_fastq_mode.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --reference_genes_gtf_file          :   Reference genes GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf).
    --params_sqanti3_qc                 :   sqanti3_qc.py parameters (default: '"--aligner_choice minimap2 --report html --force_id_ignore --isoAnnotLite"').
                                            Note that the parameters need to be wrapped in quotes.
    --params_sqanti3_filter             :   sqanti3_filter.py parameters (default: '"ml"').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header     | Description                   |
|------------|-------------------------------|
| sample_id  | Sample ID.                    |
| fastq_file | Full path to `fastq.gz` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

`--reference_genome_fasta_fai_file`
* Reference genome FASTA.FAI files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

`--reference_genes_gtf_file`
* Reference GTF files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

`--params_sqanti3_qc`
* Refer to the [Sqanti3 documentation](https://github.com/ConesaLab/SQANTI3).
* The following parameters for `sqanti3` are already included in `nexus` module for `sqanti3_qc.py` and should not be specified:
  * `--fasta`
  * `-t`
  * `-o`
  * `-d`

`--params_sqanti3_filter`
* Refer to the [Sqanti3 documentation](https://github.com/ConesaLab/SQANTI3).
* The following parameters for `sqanti3` are already included in `nexus` module for `sqanti3_filter.py` and should not be specified:
  * `--output`
  * `--dir`
