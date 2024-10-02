## novel_isoform_discovery_isoseq.nf

Discovers novel isoforms using [IsoSeq](https://isoseq.how/).

### Inputs / Outputs

| I/O    | Description                             |
|:-------|:----------------------------------------|
| Input  | Unaligned `bam` file for each sample.   | 
| Output | Clustered `bam` file for each sample.   |

### Dependencies

* `IsoSeq`

### Example

```
nexus run --nf-workflow novel_isoform_discovery_isoseq.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR
```

### Usage

```
workflow:
    1. Run isoseq 'cluster2'.

usage: nexus run --nf-workflow novel_isoform_discovery_isoseq.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'bam_file'.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`
* The `bam_file` must be an unaligned `bam` file. To generate unaligned `bam` files from `fastq` files, you can use [fastq_to_unaligned_bam.nf](/src/nexuslib/pipelines/utilities/fastq_to_unaligned_bam/) in `nexus`.
* A TSV (tab-separated values) file with the following headers:

| Header    | Description                        |
|-----------|------------------------------------|
| sample_id | Sample ID.                         |
| bam_file  | Full path to unaligned `bam` file. |

