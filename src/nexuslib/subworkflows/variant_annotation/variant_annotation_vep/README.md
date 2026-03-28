## variant_annotation_vep.nf

Annotate variants using [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html).

### Inputs / Outputs

| I/O    | Description                 |
|:-------|:----------------------------|
| Input  | `vcf` file for each sample. | 
| Output | `tsv` file for each sample. |

### Dependencies

* `VEP`

### Example

```
nexus run --nf-workflow variant_annotation_vep.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --vep_dir VEP_DIR \
    --output_dir OUTPUT_DIR \
    --params_vep '"--species homo_sapiens --database --offline --cache --assembly GRCh38"'
```

### Usage

```
workflow:
    1. Run VEP.

usage: nexus run --nf-workflow variant_annotation_vep.nf [required] [optional] [--help]

required arguments:
    -c                      :   Nextflow .config file.
    -w                      :   Nextflow work directory path.
    --samples_tsv_file      :   TSV file with the following columns: 'sample_id', 'vcf_file'.
    --vep_dir               :   VEP cache directory.
    --output_dir            :   Directory to which output files will be copied.

optional arguments:
    --params_vep            :   ClairS parameters (default: '"--species homo_sapiens --database --offline --cache --assembly GRCh38"').
                                Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header              | Description                        |
|---------------------|------------------------------------|
| sample_id           | Sample ID                          |
| tumor_bam_file      | Full path to tumor `bam` file      |
| tumor_bam_bai_file  | Full path to tumor `bam.bai` file  |
| normal_bam_file     | Full path to normal `bam` file     |
| normal_bam_bai_file | Full path to normal `bam.bai` file |

`--vep_dir`
* Path to VEP folder (e.g. `/home/users/user/.vep`). 
* Download a `tar.gz` from https://github.com/Ensembl/ensembl-vep/releases and install VEP using 
  * `perl INSTALL.pl --AUTO c --SPECIES homo_sapiens --ASSEMBLY GRCh38`
  * `perl INSTALL.pl --AUTO c --SPECIES mus_musculus --ASSEMBLY GRCm39`

`--params_vep`
* Refer to the [VEP documentation](https://useast.ensembl.org/info/docs/tools/vep/index.html).
* The following parameters for `vep` are already included in `nexus` module for `vep` and should not be specified:
  * `-i`
  * `--vcf`
  * `--dir_cache`
  * `-o`
