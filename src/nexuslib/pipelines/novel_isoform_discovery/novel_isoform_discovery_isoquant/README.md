## novel_isoform_discovery_isoquant.nf

Discovers novel isoforms using [IsoQuant](https://github.com/ablab/IsoQuant).

### Inputs / Outputs

| I/O    | Description                                                                                                     |
|:-------|:----------------------------------------------------------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                                                                                | 
| Output | IsoQuant output files including `transcript_models.gtf` and `transcript_model_reads.tsv` files for each sample. |

### Dependencies

* `IsoQuant`

### Example

```
nexus run --nf-workflow novel_isoform_discovery_isoquant.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --gtf_file GTF_FILE \
    --params_isoquant '"--data_type pacbio_ccs --sqanti_output --high_memory --complete_genedb"'
```

### Usage

```
workflow:
    1. Run isoquant.

usage: nexus run --nf-workflow novel_isoform_discovery_isoquant.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --gtf_file                      :   GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf).
    --params_isoquant               :   isoquant parameters (default: '"--data_type pacbio_ccs --sqanti_output --high_memory --complete_genedb"').
                                        Note that the parameters need to be wrapped in quotes.
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

`--gtf_file`
* Reference GTF files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

`--params_isoquant`
* Refer to the [IsoQuant documentation](https://github.com/ablab/IsoQuant).
* The following parameters for `IsoQuant` are already included in `nexus` module for `IsoQuant` and should not be specified:
  * `--reference`
  * `--genedb`
  * `--fastq_list`
  * `--threads`
  * `-o`
