## long_read_dna_variant_calling_pbsv.nf

Identifies structural DNA variants in long-read DNA BAM files using [pbsv](https://github.com/PacificBiosciences/pbsv).

### Inputs / Outputs

| I/O    | Description                  |
|:-------|:-----------------------------|
| Input  | `bam` file for each sample.  | 
| Output | `vcf` file for each sample. |

### Dependencies

* `pbsv`

### Example

```
nexus run --nf-workflow long_read_dna_variant_calling_pbsv.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --params_pbsv_discover '"--ccs --min-gap-comp-id-perc 97 --min-mapq 20"' \
    --params_pbsv_call '"--ccs --call-min-reads-per-strand-all-samples 0 --call-min-read-perc-one-sample 10 --call-min-reads-all-samples 3 --call-min-reads-one-sample 3"'
```

### Usage

```
workflow:
    1. Run Pbsv.

usage: nexus run --nf-workflow long_read_dna_variant_calling_pbsv.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --params_pbsv_discover          :   Pbsv 'discover' parameters (default: '"--ccs --min-gap-comp-id-perc 97 --min-mapq 20"').
                                        Note that the parameters need to be wrapped in quotes.
    --params_pbsv_call              :   Pbsv 'call' parameters (default:
                                        '"--ccs
                                          --call-min-reads-per-strand-all-samples 0
                                          --call-min-read-perc-one-sample 10
                                          --call-min-reads-all-samples 3
                                          --call-min-reads-one-sample 3"').
                                        Note that the parameters need to be wrapped in quotes.
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`

| Header       | Description                 |
|--------------|-----------------------------|
| sample_id    | Sample ID                   |
| bam_file     | Full path to `bam` file     |
| bam_bai_file | Full path to `bam.bai` file |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--params_pbsv_discover`
* Refer to the [pbsv documentation](https://github.com/PacificBiosciences/pbsv).

`--params_pbsv_call`
* Refer to the [pbsv documentation](https://github.com/PacificBiosciences/pbsv).
* The following parameters for `pbsv` are already included in `nexus` module for `pbsv call` and should not be specified:
  * `--num-threads`
