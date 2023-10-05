## pacbio_dna_structural_variants.nf

Identifies structural DNA variants with long-read DNA sequencing 
(BAM) files using 
[Sniffles2](https://github.com/fritzsedlazeck/Sniffles), 
[Pbsv](https://github.com/PacificBiosciences/pbsv), 
[SVIM](https://github.com/eldariont/svim), 
[cuteSV](https://github.com/tjiangHIT/cuteSV).

### Inputs / Outputs

| Input(s)                  | Output(s)  |
|---------------------------|------------|
| `BAM` files | `VCF` files |

### Dependencies

* `sniffles2`
* `pbsv`
* `svim`
* `cutesv`

### Usage

```shell
nextflow run pacbio_dna_structural_variants.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                    :   Directory to which output files will be copied.
    --reference_genome_fasta_file   :   Reference genome FASTA file.
    --sniffles2                     :   Sniffles2 path.
    --sniffles2_params              :   Sniffles2 parameters (e.g. "").
    --pbsv                          :   pbsv path.
    --pbsv_discover_params          :   pbsv 'discover' parameters (e.g. "").
    --pbsv_call_params              :   pbsv 'call' parameters (e.g. "").
    --svim                          :   SVIM path.
    --svim_params                   :   SVIM parameters (e.g. "").
    --cutesv                        :   cuteSV path.
    --cutesv_params                 :   cuteSV parameters (e.g. "").

optional arguments:
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`--sample_tsv_file`

| Header       | Description                 |
|--------------|-----------------------------|
| sample_id    | Sample ID                   |
| bam_file     | Full path to `BAM` file     |
| bam_bai_file | Full path to `BAM.BAI` file |
