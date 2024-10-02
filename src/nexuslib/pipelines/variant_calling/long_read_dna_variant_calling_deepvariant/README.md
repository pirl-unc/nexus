## long_read_dna_variant_calling_deepvariant.nf

Identifies small DNA variants (SNVs and INDELs) in long-read DNA BAM files using [DeepVariant](https://github.com/google/deepvariant/).

### Inputs / Outputs

| I/O    | Description                                |
|:-------|:-------------------------------------------|
| Input  | `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.                |

### Dependencies

* `DeepVariant`

### Example

```
nexus run --nf-workflow long_read_dna_variant_calling_deepvariant.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --containerization singularity \
    --deepvariant_model_type PACBIO \
    --deepvariant_bin_path /opt/deepvariant/bin/run_deepvariant\
    --deepvariant_bin_version 1.6.0 \
    --deepvariant_input_path INPUT_PATH \
    --deepvariant_output_path OUTPUT_PATH
```

### Usage

```
workflow:
    1. Run DeepVariant.

usage: nexus run --nf-workflow long_read_dna_small_variants.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns: 
                                            'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --containerization                  :   Containerization ('singularity' or 'docker'; default: 'singularity').
    --deepvariant_model_type            :   DeepVariant --model_type parameter value (default: 'PACBIO').
    --deepvariant_bin_path              :   DeepVariant bin path (default: '/opt/deepvariant/bin/run_deepvariant').
    --deepvariant_bin_version           :   DeepVariant bin version (default: '1.6.0').
    --deepvariant_input_path            :   DeepVariant input path (default: '/datastore/').
    --deepvariant_output_path           :   DeepVariant output path (default: '/datastore/').
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                    |
|--------------|--------------------------------|
| sample_id    | Sample ID.                     |
| bam_file     | Full path to `bam` file.       |
| bam_bai_file | Full path to `bam.bai` file.   |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

