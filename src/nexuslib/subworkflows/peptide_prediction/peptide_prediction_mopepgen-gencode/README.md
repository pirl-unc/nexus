## peptide_prediction_mopepgen-gencode.nf

Discovers mutant peptide sequences using [mopepgen](https://github.com/uclahs-cds/package-moPepGen).

### Inputs / Outputs

| I/O    | Description                                                                                                                  |
|:-------|:-----------------------------------------------------------------------------------------------------------------------------|
| Input  | `mutect2`, `strelka2`, `delly2`, `manta`, `dysgu`, `reditools2`, `arriba`, `rmats`, `circexplorer2` outputs for each sample. | 
| Output | moPepGen outputs for each sample.                                                                                            |

### Dependencies

* `mopepgen`

### Example

```
nexus run --nf-workflow peptide_prediction_mopepgen-gencode.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --reference_proteome_fasta_file REFERENCE_PROTEOME_FASTA_FILE \
    --vep_dir VEP_DIR \
    --params_vep '"--species homo_sapiens --offline --cache --assembly GRCh38"' \
    --params_filter_vep '"--filter Source = GENCODE"' \
    --params_mopepgen_parsevep '""' \
    --params_mopepgen_parsereditools '""' \
    --params_mopepgen_parsearriba '""' \
    --params_mopepgen_parsermats '""' \
    --params_mopepgen_parsecircexplorer2 '""' \
    --params_mopepgen_callvariant '""'    
```

### Usage

```
workflow:
    1. Run moPepGen callVariant.

usage: nexus run --nf-workflow peptide_prediction_mopepgen-gencode.nf [required] [optional] [--help]

required arguments:
    -c                                      :   Nextflow .config file.
    -w                                      :   Nextflow work directory path.
    --samples_tsv_file                      :   TSV file with the following columns: 'sample_id',
                                                'mutect2_vcf_file',
                                                'strelka2_snv_vcf_file',
                                                'strelka2_indel_vcf_file',
                                                'delly2_vcf_file',
                                                'manta_vcf_file',
                                                'dysgu_vcf_file',
                                                'reditools2_tsv_file',
                                                'arriba_tsv_file',
                                                'rmats_output_dir',
                                                'circexplorer2_tsv_file'.
                                                Use 'NA' or empty string for missing files.
    --output_dir                            :   Directory to which output files will be copied.
    --reference_genome_fasta_file           :   Reference genome FASTA file (GENCODE).
    --reference_genes_gtf_file              :   Reference genes GTF file (GENCODE).
    --reference_proteome_fasta_file         :   Reference proteome FASTA file (GENCODE).
    --vep_dir                               :   VEP cache directory.

optional arguments:
    --params_vep                            :   VEP parameters (default: "--species homo_sapiens --offline --cache --assembly GRCh38").
                                                Note that the parameters need to be wrapped in quotes.
    --params_filter_vep                     :   filter_vep parameters (default: '"--filter Source = GENCODE"').
                                                Note that the parameters need to be wrapped in quotes.
    --params_mopepgen_parsevep              :   moPepGen parseVEP parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
    --params_mopepgen_parsereditools        :   moPepGen parseREDItools parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
    --params_mopepgen_parsearriba           :   moPepGen parseArriba parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
    --params_mopepgen_parsermats            :   moPepGen parseRMATS parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
    --params_mopepgen_parsecircexplorer2    :   moPepGen parseCIRCexplorer parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
    --params_mopepgen_callvariant           :   moPepGen callVariant parameters (default: '""').
                                                Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header                   | Description                                                                           |
|--------------------------|---------------------------------------------------------------------------------------|
| sample_id                | Sample ID.                                                                            |
| mutect2_vcf_file         | Full path to Mutect2 'vcf` file. Use 'NA' or empty string "" for missing file.        |
| strelka2_snv_vcf_file    | Full path to Strelka2 SNV 'vcf` file. Use 'NA' or empty string "" for missing file.   |
| strelka2_indel_vcf_file  | Full path to Strelka2 INDEL 'vcf` file. Use 'NA' or empty string "" for missing file. |
| delly2_vcf_file          | Full path to Delly2 'vcf` file. Use 'NA' or empty string "" for missing file.         |
| manta_vcf_file           | Full path to Manta 'vcf` file. Use 'NA' or empty string "" for missing file.          |
| dysgu_vcf_file           | Full path to Dysgu 'vcf` file. Use 'NA' or empty string "" for missing file.          |
| reditools2_tsv_file      | Full path to REDItool2 'tsv` file. Use 'NA' or empty string "" for missing file.      |
| arriba_tsv_file          | Full path to Arriba 'tsv` file. Use 'NA' or empty string "" for missing file.         |
| rmats_output_dir         | Full path to RMATS output directory. Use 'NA' or empty string "" for missing file.    |
| circexplorer2_tsv_file   | Full path to CIRCexplorer2 'tsv` file. Use 'NA' or empty string "" for missing file.  |

`--reference_genome_fasta_file`
* GENCODE reference genome FASTA (`.fasta`) file. An uncompressed `.fasta` file should be supplied.

`--reference_genes_gtf_file`
* GENCODE reference genes GTF (`.gtf`) file. An uncompressed `.gtf` file should be supplied.

`--reference_proteome_fasta_file`
* GENCODE reference proteome FASTA (`.fasta`) file. An uncompressed `.fasta` file should be supplied.

`--vep_dir`
* VEP cache directory. Download a `tar.gz` from https://github.com/Ensembl/ensembl-vep/releases and install VEP using `INSTALL.pl`. This should create a VEP cache directory.

`--params_vep`
* VEP parameters (default: `'"--species homo_sapiens --offline --cache --assembly GRCh38 --distance 0"'`).

`--params_filter_vep`
* filter_vep parameters (default: `'"--filter Source = GENCODE"'`).

`--params_mopepgen_parsevep`
* moPepGen parseVEP parameters (default: `'""'`)
* Refer to the [mopepgen documentation](https://uclahs-cds.github.io/package-moPepGen/parse-vep/).
* The following parameters for `moPepGen parseVEP` are already included in `nexus` module for `mopepgen` and should not be specified:
  * `-i`
  * `--output-path`
  * `--genome-fasta`
  * `--annotation-gtf`
  * `--reference-source`
  * `--source`

`--params_mopepgen_parsereditools`
* moPepGen parseREDItools parameters (default: `'""'`)
* Refer to the [mopepgen documentation](https://uclahs-cds.github.io/package-moPepGen/parse-vep/).
* The following parameters for `moPepGen parseREDItools` are already included in `nexus` module for `mopepgen` and should not be specified:
  * `-i`
  * `--output-path`
  * `--annotation-gtf`
  * `--reference-source`
  * `--source`

`--params_mopepgen_parsearriba`
* moPepGen parseArriba parameters (default: `'""'`)
* Refer to the [mopepgen documentation](https://uclahs-cds.github.io/package-moPepGen/parse-vep/).
* The following parameters for `moPepGen parseArriba` are already included in `nexus` module for `mopepgen` and should not be specified:
  * `-i`
  * `--output-path`
  * `--genome-fasta`
  * `--annotation-gtf`
  * `--reference-source`
  * `--source arriba`

`--params_mopepgen_parsermats`
* moPepGen parseRMATS parameters (default: `'""'`)
* Refer to the [mopepgen documentation](https://uclahs-cds.github.io/package-moPepGen/parse-vep/).
* The following parameters for `moPepGen parseArriba` are already included in `nexus` module for `mopepgen` and should not be specified:
  * `--se`
  * `--a5ss`
  * `--a3ss`
  * `--mxe`
  * `--ri`
  * `--output-path`
  * `--genome-fasta`
  * `--annotation-gtf`
  * `--reference-source`
  * `--source rmats`

`--params_mopepgen_parsecircexplorer2`
* moPepGen parseCIRCexplorer parameters (default: `'""'`)
* Refer to the [mopepgen documentation](https://uclahs-cds.github.io/package-moPepGen/parse-vep/).
* The following parameters for `moPepGen parseCIRCexplorer` are already included in `nexus` module for `mopepgen` and should not be specified:
  * `-i`
  * `--output-path`
  * `--annotation-gtf`
  * `--reference-source`
  * `--source circexplorer2`

`--params_mopepgen_callvariant`
* moPepGen callVariant parameters (default: `'""'`)
* Refer to the [mopepgen documentation](https://uclahs-cds.github.io/package-moPepGen/parse-vep/).
* The following parameters for `moPepGen callVariant` are already included in `nexus` module for `mopepgen` and should not be specified:
  * `-i`
  * `--output-path`
  * `--threads`
  * `--genome-fasta`
  * `--annotation-gtf`
  * `--proteome-fasta`
