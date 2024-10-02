## antigen_presentation_prediction_mhcflurry2.nf

Predict antigen presentation using [MHCflurry 2.0](https://github.com/openvax/mhcflurry).

### Inputs / Outputs

| I/O    | Description                                |
|:-------|:-------------------------------------------|
| Input  | MHCflurry 2.0 `csv` file for each sample.  | 
| Output | MHCflurry 2.0 `csv` file for each sample.  |

### Dependencies

* `MHCflurry 2.0`

### Example
```
nexus run --nf-workflow antigen_presentation_prediction_mhcflurry2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --params_mhcflurry2_predict '""' \
```

### Usage

```
workflow:
    1. Run mhcflurry-predict command.

usage: nexus run --nf-workflow antigen_presentation_prediction_mhcflurry2.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id',
                                            'mhcflurry2_input_csv_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --params_mhcflurry2_predict         :   mhcflurry-predict parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow `config` file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header                    | Description                               |
|---------------------------|-------------------------------------------|
| sample_id                 | Sample ID.                                |
| mhcflurry2_input_csv_file | Full path to MHCflurry2 input `csv` file. |

`--params_mhcflurry2_predict`
* Refer to the [MHCflurry2 documentation](https://github.com/openvax/mhcflurry).
* The following parameters for `mhcflurry-predict` are already included in `nexus` module and should not be specified:
  * `--out`
