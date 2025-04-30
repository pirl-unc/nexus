## novel_isoform_discovery_isotools.nf

Discovers novel isoforms using [IsoTools](https://isotools.readthedocs.io/en/latest/).

### Inputs / Outputs

| I/O    | Description                                |
|:-------|:-------------------------------------------|
| Input  | `bam` and `bam.bai` files for each sample. | 
| Output | IsoTools output files for each sample.     |

### Dependencies

* `IsoTools`

### Example

```
nexus run --nf-workflow novel_isoform_discovery_isotools.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --params_isotools '""'
```

### Usage

```
workflow:
    1. Run Isotools.

usage: nexus run --nf-workflow novel_isoform_discovery_isotools.nf [required] [optional] [--help]

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
    --params_isotools                   :   run_isotools parameters (default: '""').
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

`--params_isotools`
* Refer to the [isotools documentation](https://isotools.readthedocs.io/en/latest/).
* The following parameters for `isotools` are already included in `nexus` module for `run_isotools` and should not be specified:
  * `--anno`
  * `--genome`
  * `--samples`
  * `--file_prefix`
  * `--gtf_out`
