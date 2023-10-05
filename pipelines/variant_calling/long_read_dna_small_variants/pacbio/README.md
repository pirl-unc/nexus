## pacbio_dna_small_variants.nf

Identifies small DNA variants (SNVs and INDELs) with long-read DNA sequencing 
(BAM) files using [DeepVariant](https://github.com/google/deepvariant/).

### Inputs / Outputs

| Input(s)                  | Output(s)  |
|---------------------------|------------|
| `BAM` files | `VCF` files |

### Dependencies

* `deepvariant`

### Usage

```shell
nextflow run pacbio_dna_small_variants.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'bam_file', 'bam_bai_file'
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --singularity                       :   Singularity path.
    --deepvariant_bin_path              :   DeepVariant bin path.
    --deepvariant_bin_version           :   DeepVariant bin version.
    --deepvariant_lib_path              :   DeepVariant lib path.
    --deepvariant_model_type            :   DeepVariant --model_type parameter value.

optional arguments:
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`--sample_tsv_file`

| Header       | Description                 |
|--------------|-----------------------------|
| sample_id    | Sample ID                   |
| bam_file     | Full path to `BAM` file     |
| bam_bai_file | Full path to `BAM.BAI` file |
