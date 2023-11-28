## paired_end_read_rna_quantification_salmon_mapping_mode.nf

Quantifies RNA abundance with paired-end RNA sequencing (FASTQ) files using 
[salmon](https://combine-lab.github.io/salmon/) in mapping mode.

### Inputs / Outputs

| Input(s)                    | Output(s)                                                                                         |
| --------------------------- |---------------------------------------------------------------------------------------------------|
| `R1 and R2 FASTQ.GZ` files  | `quant.sf` files |

### Dependencies

* `salmon`

### Usage

```shell
nextflow run paired_end_read_rna_quantification_salmon_mapping_mode.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file_1', 'fastq_file_2', 'fasta_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --salmon                            :   salmon path.
    --salmon_index_params               :   salmon 'index' parameters (e.g. '--gencode').
    --salmon_quant_params               :   salmon 'quant' parameters
                                            (e.g. '--libType <library_type> --geneMap <gtf_file> --seqBias --gcBias --posBias').
    --salmon_index                      :   salmon index path.

optional arguments:
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`--sample_tsv_file`

| Header       | Description                     |
| ------------ | ------------------------------- |
| sample_id    | Sample ID                       |
| fastq_file_1 | Full path to R1 `FASTQ.GZ` file |
| fastq_file_2 | Full path to R2 `FASTQ.GZ` file |
| fasta_file   | Full path to transcripts `FASTA` file |

`--salmon_quant_params`

Refer to the [salmon documentation](https://combine-lab.github.io/salmon/).

`--salmon_index`

A [salmon index](https://combine-lab.github.io/salmon/) needs to be generated prior to running this workflow script.
