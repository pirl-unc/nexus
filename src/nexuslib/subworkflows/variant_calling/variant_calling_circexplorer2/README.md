## variant_calling_circexplorer2.nf

Identifies circular RNAs in paired-end RNA-seq files using [CIRCexplorer2](https://circexplorer2.readthedocs.io/en/latest/).

### Inputs / Outputs

| I/O    | Description                               |
|:-------|:------------------------------------------|
| Input  | CIRCexplorer2 input file for each sample. | 
| Output | `VCF` file for each sample.             |

### Dependencies

* `CIRCexplorer2`

### Example

```
nexus run --nf-workflow variant_calling_circexplorer2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --circexplorer2_gene_annotation_txt_file CIRCEXPLORER2_GENE_ANNOTATION_TXT_FILE \
    --params_circexplorer2_parse '"-t STAR"' \
    --params_circexplorer2_annotate '""'
```

### Usage

```
workflow:
    1. Run CIRCexplorer2 parse.
    2. Run CIRCexplorer2 annotate.

usage: nexus run --nf-workflow variant_calling_circexplorer2.nf [required] [optional] [--help]

required arguments:
    -c                                          :   Nextflow .config file.
    -w                                          :   Nextflow work directory path.
    --samples_tsv_file                          :   TSV file with the following columns:
                                                    'sample_id', 'input_file'.
    --output_dir                                :   Directory to which output files will be copied.
    --reference_genome_fasta_file               :   Reference genome FASTA file.
    --circexplorer2_gene_annotation_txt_file    :   CIRCexplorer2 gene annotation TXT file.

optional arguments:
    --params_circexplorer2_parse                :   CIRCexplorer2 parse parameters (default: '"-t STAR"').
                                                    Note that the parameters need to be wrapped in quotes.
    --params_circexplorer2_annotate             :   CIRCexplorer2 annotate parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header     | Description                                                    |
|------------|----------------------------------------------------------------|
| sample_id  | Sample ID                                                      |
| input_file | Full path to input file (e.g. STAR Chimeric.out.junction file) |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta` or `.fasta.gz`) file.

`--circexplorer2_gene_annotation_txt_file`
* CIRCexplorer gene annotation TXT file. This file can be downloaded by running `fetch_ucsc.py hg38 kg hg38_ens.txt`. `fetch_ucsc.py` is available as part of CIRCexplorer2.

`--params_circexplorer2_parse`
* Refer to the [CIRCexplorer2 documentation](https://circexplorer2.readthedocs.io/en/latest/).
* The following parameters for `CIRCexplorer2 parse` are already included in `nexus` module for `CIRCexplorer2` and should not be specified:
  * `-b`

`--params_circexplorer2_annotate`
* Refer to the [CIRCexplorer2 documentation](https://circexplorer2.readthedocs.io/en/latest/).
* The following parameters for `CIRCexplorer2 annotate` are already included in `nexus` module for `CIRCexplorer2` and should not be specified:
  * `-r`
  * `-g`
  * `-b`
  * `-o`
