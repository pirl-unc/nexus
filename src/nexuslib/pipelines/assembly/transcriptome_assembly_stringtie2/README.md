## transcriptome_assembly_stringtie2.nf

Assemble transcript sequences using long RNA reads using [StringTie2](https://github.com/skovaka/stringtie2).

### Inputs / Outputs

| I/O    | Description                      |
|:-------|:---------------------------------|
| Input  | `fastq.gz` file for each sample. | 
| Output | `gtf` file for each sample.      |

### Dependencies

* `RNA-Bloom2`

### Example

```
nexus run --nf-workflow transcriptome_assembly_stringtie2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --params_stringtie2 '"-L"'
```

### Usage

```
workflow:
    1. Assemble transcripts using StringTie2.

usage: nexus run --nf-workflow transcriptome_assembly_stringtie2.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --gtf_file                          :   GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf.gz).
    --params_stringtie2                 :   StringTie2 parameters (default: '"-L"').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                     |
|--------------|---------------------------------|
| sample_id    | Sample ID.                      |
| fastq_file   | Full path to `fastq.gz` file    |

`--params_stringtie2`
* Refer to the [StringTie2 documentation](https://github.com/skovaka/stringtie2).
* The following parameters for `StringTie2` are already included in `nexus` module for `StringTie2` and should not be specified:
  * `-o`
  * `-G`
  * `-p`
