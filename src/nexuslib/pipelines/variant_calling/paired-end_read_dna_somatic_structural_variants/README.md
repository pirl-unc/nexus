## paired_end_human_dna_somatic_small_variants.nf

Identifies structural DNA variants with paired-end DNA sequencing 
(BAM) files using 
[Delly2](https://github.com/dellytools/delly), 
[Lumpy](https://github.com/arq5x/lumpy-sv).

### Inputs / Outputs

| Input(s)                  | Output(s)  |
|---------------------------|------------|
| `BAM` files | `VCF` files |

### Dependencies

* `delly`
* `lumpyexpress`

### Usage

```shell
nextflow run paired_end_dna_somatic_structural_variants.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file                          :   TSV file with the following columns:
                                                    'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file', 'tumor_sample_id', 'normal_sample_id'
    --output_dir                                :   Directory to which output files will be symlinked.
    --reference_genome_fasta_file               :   Reference genome FASTA file.
    --delly2                                    :   DELLY2 path.
    --delly2_call_params                        :   DELLY2 'call' parameters (e.g. "").
    --bcftools                                  :   bcftools path.
    --python2                                   :   python2 path.
    --lumpy_express                             :   LUMPY Express path.
    --lumpy_extract_split_reads_script_file     :   LUMPY 'extractSplitReads_BwaMem' file.
    --lumpy_config_file                         :   LUMPY config file.
    --samtools                                  :   samtools path.
```

### Parameters

`--sample_tsv_file`

| Header       | Description                 |
|--------------|-----------------------------|
| sample_id    | Sample ID                   |
| bam_file     | Full path to `BAM` file     |
| bam_bai_file | Full path to `BAM.BAI` file |
