## novel_isoform_discovery_flair.nf

Discovers novel isoforms using [flair](https://github.com/BrooksLabUCSC/flair).

### Inputs / Outputs

| I/O    | Description                                                         |
|:-------|:--------------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                                    | 
| Output | `fasta`, `gtf`, `bed`, `bam`, and `bam.bai` files  for each sample. |

### Dependencies

* `flair`

### Usage

```
workflow:
    1. Run flair align.
    2. Run flair correct.
    3. Run flair collapse.

usage: nexus run --nf-workflow novel_isoform_discovery_flair.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --gtf_file                      :   Reference transcriptome GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf).
    --flair                         :   flair path (default: flair).
    --flair_align_params            :   flair 'align' parameters (default: '" "').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
    --flair_correct_params          :   flair 'correct' parameters (default: '" "').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
    --flair_collapse_params         :   flair 'collapse' parameters (default: '" "').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/configs)

`--sample_tsv_file`

| Header     | Description                   |
| ---------- |-------------------------------|
| sample_id  | Sample ID.                    |
| fastq_file | Full path to `fastq.gz` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.
