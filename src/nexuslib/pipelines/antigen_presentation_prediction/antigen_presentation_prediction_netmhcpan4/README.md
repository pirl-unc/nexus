## antigen_presentation_prediction_netmhcpan4.nf

Predict antigen presentation using [NetMHCpan4.x](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/).

### Inputs / Outputs

| I/O    | Description                                   |
|:-------|:----------------------------------------------|
| Input  | `fasta` file and MHC alleles for each sample. | 
| Output | NetMHCpan 4.x `txt` file for each sample.     |

### Dependencies

* `NetMHCpan4.x`

### Example
```
nexus run --nf-workflow antigen_presentation_prediction_netmhcpan4.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --netmhcpan_home_dir NETMHCPAN_HOME_DIR \
    --netmhcpan NETMHCPAN \
    --params_netmhcpan '"-BA -s -l 8,9,10,11"' \
```

### Usage

```
workflow:
    1. Run netMHCpan command.

usage: nexus run --nf-workflow antigen_presentation_prediction_netmhcpan4.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --netmhcpan_home_dir                :   NetMHCpan home directory path.
    --netmhcpan                         :   NetMHCpan executable path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id',
                                            'fasta_file',
                                            'mhc_alleles'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --params_netmhcpan                  :   netmhcpan parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow `config` file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header     | Description               |
|------------|---------------------------|
| sample_id  | Sample ID.                |
| fasta_file | Full path to `fasta` file. |
| mhc_alleles | MHC alleles to predict.   |

`--params_netmhcpan`
* Refer to the [NetMHCpan4 documentation](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/).
* The following parameters for `netMHCpan` are already included in `nexus` module and should not be specified:
  * `-a`
  * `-f`
  * `-inptype`

### Notes

After you run this workflow, you can use a `Nexus` built-in executable `convert_netmhcpan_txt_to_tsv` 
(installed with `Nexus`) to convert the `NetMHCpan4.x` output `txt` file to a `tsv` file:

```
convert_netmhcpan_txt_to_tsv [-h] 
  [-i INPUT_TXT_FILE]
  [-o OUTPUT_TSV_FILE]
```
