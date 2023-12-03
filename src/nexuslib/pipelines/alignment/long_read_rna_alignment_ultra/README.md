## long_read_rna_alignment_ultra.nf

Aligns long RNA reads (Pacific Biosciences or Oxford Nanopore) to a reference genome using [uLTRA](https://github.com/ksahlin/ultra).

### Inputs / Outputs

| I/O    | Description                                                     |
|:-------|:----------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                                | 
| Output | MD-tagged and sorted `bam` and `bam.bai` files for each sample. |

### Dependencies

* `uLTRA`
* `samtools`

### Usage

```
workflow:
    1. Align reads (fastq.gz files) to a reference genome using uLTRA.
    2. Generate MD tags.
    3. Sort MD-tagged bam file.

usage: nexus run --nf-workflow long_read_rna_alignment_ultra.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'fastq_file'.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --ultra                         :   uLTRA path (default: uLTRA).
    --ultra_index                   :   uLTRA index path (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/ultra/hg38_index/).
    --ultra_params                  :   uLTRA parameters (default: '"--isoseq "').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
    --samtools                      :   samtools path (default: samtools).
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/configs)

`--sample_tsv_file`

| Header     | Description                  |
| ---------- |------------------------------|
| sample_id  | Sample ID.                   |
| fastq_file | Full path to `fastq.gz` file |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on the LBG.

`--ultra_index`
* An [uLTRA index](https://github.com/ksahlin/ultra) needs to be generated prior to running this workflow. 
* Prebuilt uLTRA indices are available in `/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/ultra/` on the LBG.

`--ultra_params`
* Refer to the [uLTRA documentation](https://github.com/ksahlin/ultra).