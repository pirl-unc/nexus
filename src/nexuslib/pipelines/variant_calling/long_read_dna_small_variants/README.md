## long_read_dna_small_variants.nf

Identifies small DNA variants (SNVs and INDELs) with long-read DNA sequencing 
(BAM) files using [DeepVariant](https://github.com/google/deepvariant/).

### Inputs / Outputs

| I/O    | Description                   |
|:-------|:------------------------------|
| Input  | `bam` file for each sample.   | 
| Output | `vcf` file for each sample.   |

### Dependencies

* `deepvariant`

### Usage

```
workflow:
    1. Run DeepVariant.

usage: nexus run --nf-workflow long_read_dna_small_variants.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --containerization                  :   Containerization ('singularity' or 'docker'; default: 'singularity').
    --deepvariant_bin_path              :   DeepVariant bin path (e.g. '/opt/deepvariant/bin/run_deepvariant').
    --deepvariant_bin_version           :   DeepVariant bin version (e.g. 1.6.0).
    --deepvariant_input_path            :   DeepVariant input path (e.g. /path/to/input/).
    --deepvariant_output_path           :   DeepVariant output path (e.g. /path/to/output/).
    --deepvariant_model_type            :   DeepVariant --model_type parameter value (e.g. 'WGS').

optional arguments:
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/configs)

`--sample_tsv_file`

| Header       | Description                    |
|--------------|--------------------------------|
| sample_id    | Sample ID.                     |
| bam_file     | Full path to `bam` file.       |
| bam_bai_file | Full path to `bam.bai` file.   |

`--reference_genome_fasta_file`

* Reference genome FASTA files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/references/` on the LBG.
