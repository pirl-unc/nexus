## long_read_rna_quantification_liqa.nf

Quantify RNA abundance in long RNA reads using [LIQA](https://github.com/WGLab/LIQA).

### Inputs / Outputs

| I/O    | Description                                |
|:-------|:-------------------------------------------|
| Input  | `bam` and `bam.bai` files for each sample. | 
| Output | LIQA output files                          |

### Dependencies

* `LIQA`

### Example

```
nexus run --nf-workflow long_read_rna_quantification_liqa.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --liqa_refgene_file LIQA_REFGENE_FILE \
    --output_dir OUTPUT_DIR \
    --params_liqa_quantify PARAMS_LIQA_QUANTIFY
```

### Usage

```
workflow:
    1. Quantify RNA in long-read FASTQ files using LIQA.

usage: nexus run --nf-workflow long_read_rna_quantification_liqa.nf [required] [optional] [--help]

required arguments:
    -c                                      :   Nextflow .config file.
    -w                                      :   Nextflow work directory path.
    --samples_tsv_file                      :   TSV file with the following columns:
                                                'sample_id', 'bam_file', 'bam_bai_file'.
    --liqa_refgene_file                     :   LIQA refgene file.
    --output_dir                            :   Directory to which output files will be copied.

optional arguments:
    --params_liqa_quantify                  :   LIQA quantify parameters (default: '"-max_distance 20 -f_weight 1"').
                                                Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                       :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--liqa_refgene_file`
* LIQA refgene files can be found in 
`/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/liqa/` on LBG.

`--params_liqa_quantify`
* Refer to the [LIQA documentation](https://github.com/WGLab/LIQA/blob/master/doc/Usage.md).
* The following parameters for `LIQA` are already included in `nexus` module for `liqa -task quantify` and should not be specified:
  * `-refgene`
  * `-bam`
  * `-out`

