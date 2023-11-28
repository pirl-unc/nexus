## pacbio_dna_somatic_structural_variants.nf

Identifies somatic structural DNA variants with long-read DNA sequencing 
(BAM) files using 
[Savana](https://github.com/cortes-ciriano-lab/savana).

### Inputs / Outputs

| Input(s)                  | Output(s)  |
|---------------------------|------------|
| `BAM` files | `VCF` files |

### Dependencies

* `savana`

### Usage

```shell
nextflow run pacbio_dna_somatic_structural_variants.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'.
    --output_dir                    :   Directory to which output files will be copied.
    --reference_genome_fasta_file   :   Reference genome FASTA file.
    --savana                        :   Savana path.
    --savana_params                 :   Savana parameters (e.g. "").

optional arguments:
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`--sample_tsv_file`

| Header              | Description                        |
|---------------------|------------------------------------|
| sample_id           | Sample ID                          |
| tumor_bam_file      | Full path to tumor `BAM` file      |
| tumor_bam_bai_file  | Full path to tumor `BAM.BAI` file  |
| normal_bam_file     | Full path to normal `BAM` file     |
| normal_bam_bai_file | Full path to normal `BAM.BAI` file |
