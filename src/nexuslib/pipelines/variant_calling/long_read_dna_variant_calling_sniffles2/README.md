## long_read_dna_variant_calling_sniffles2.nf

Identifies structural DNA variants in long-read DNA BAM files using [Sniffles2](https://github.com/fritzsedlazeck/Sniffles).

### Inputs / Outputs

| I/O    | Description                  |
|:-------|:-----------------------------|
| Input  | `bam` file for each sample.  | 
| Output | `vcf` file for each sample. |

### Dependencies

* `cuteSV`

### Example

```
nexus run --nf-workflow long_read_dna_variant_calling_sniffles2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --params_sniffles2 '"--minsupport 3 --minsvlen 30 --mapq 20 --output-rnames"'
```

### Usage

```
workflow:
    1. Run Sniffles2.

usage: nexus run --nf-workflow long_read_dna_variant_calling_sniffles2.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --params_sniffles2                  :   Sniffles2 parameters (default: '"--minsupport 3 --minsvlen 30 --mapq 20 --output-rnames"').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
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

`--reference_genome_fasta_fai_file`
* Reference genome FASTA.FAI files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--params_sniffles2`
* Refer to the [Sniffles2 documentation](https://github.com/fritzsedlazeck/Sniffles).
* The following parameters for `Sniffles2` are already included in `nexus` module for `sniffles` and should not be specified:
  * `--input`
  * `--vcf`
  * `--snf`
  * `--sample-id`
  * `--reference`
  * `--threads`