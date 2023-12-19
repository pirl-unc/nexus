## novel_isoform_discovery_isoquant.nf

Discovers novel isoforms using [isoquant](https://github.com/ablab/IsoQuant).

### Inputs / Outputs

| I/O    | Description                                                                                     |
|:-------|:------------------------------------------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                                                                | 
| Output | Various `isoquant` output files including `transcript_models.gtf`, `transcript_model_reads.tsv`. |

### Dependencies

* `isoquant`

### Usage

```
workflow:
    1. Run isoquant.

usage: nexus run --nf-workflow novel_isoform_discovery_isoquant.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --isoquant                      :   isoquant path (default: isoquant.py).
    --isoquant_params               :   isoquant parameters (default: '"--data_type pacbio_ccs --sqanti_output --high_memory --complete_genedb --genedb /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf "').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/configs)

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on LBG.

`--isoquant_params`
* Make sure to include `--fastq_list experiments.txt`.
For the format of this file, refer to [here](https://github.com/ablab/IsoQuant).
