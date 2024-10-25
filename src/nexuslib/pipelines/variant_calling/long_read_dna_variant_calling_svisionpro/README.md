## long_read_dna_variant_calling_svisionpro.nf

Identifies structural DNA variants in long-read DNA BAM files using [SVision-pro](https://github.com/songbowang125/SVision-pro).

### Inputs / Outputs

| I/O    | Description                                                |
|:-------|:-----------------------------------------------------------|
| Input  | Tumor and normal `bam` and `bam.bai` files for each sample. | 
| Output | `vcf` file for each sample.                                |

### Dependencies

* `SVision-pro`

### Example

```
nexus run --nf-workflow long_read_dna_variant_calling_svisionpro.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --svisionpro_model_file SVISIONPRO_MODEL_FILE \
    --params_svisionpro '"--detect_mode somatic --preset hifi --min_supp 3 --min_mapq 20 --min_sv_size 30 --max_sv_size 1000000 --device cpu --img_size 256"' \
    --params_svisionpro_extract '"--extract somatic --min_supp 3"'
```

### Usage

```
workflow:
    1. Run SVision-pro.

usage: nexus run --nf-workflow long_read_dna_variant_calling_svisionpro.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --svisionpro_model_file             :   SVision-pro model file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/svisionpro/SVision-pro-1.8/src/pre_process/model_liteunet_256_8_16_32_32_32.pth).
    --params_svisionpro                 :   SVision-pro parameters (default: '"--detect_mode somatic --preset hifi --min_supp 3 --min_mapq 20 --min_sv_size 30 --max_sv_size 1000000 --device cpu --img_size 256"').
                                            Note that the parameters need to be wrapped in quotes.
    --params_svisionpro_extract         :   SVision-pro extract_op.py parameters (default: '"--extract somatic --min_supp 3"').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
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

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--reference_genome_fasta_fai_file`
* Reference genome FASTA.FAI files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--svisionpro_model_file`
* Refer to the [SVision-pro documentation](https://github.com/songbowang125/SVision-pro).

`--params_svisionpro`
* Refer to the [SVision-pro documentation](https://github.com/songbowang125/SVision-pro).
* The following parameters for `SVision-pro` are already included in `nexus` module for `SVision-pro` and should not be specified:
  * `--target_path`
  * `--base_path`
  * `--genome_path`
  * `--model_path`
  * `--out_path`
  * `--sample_name`

`--params_svisionpro_extract`
* Refer to the [SVision-pro documentation](https://github.com/songbowang125/SVision-pro).
* The following parameters for `extract_op.py` are already included in `nexus` module for `SVision-pro` and should not be specified:
  * `--input_vcf`
