## paired-end_read_dna_variant_calling_svaba.nf

Identifies somatic structural DNA variants in paired-read DNA BAM files using [Svaba](https://github.com/walaj/svaba).

### Inputs / Outputs

| I/O    | Description                                   |
|:-------|:----------------------------------------------|
| Input  | Tumor and normal `bam` files for each sample. | 
| Output | `vcf` file for each sample.                   |

### Dependencies

* `Svaba`

### Example

```
nexus run --nf-workflow paired-end_read_dna_variant_calling_svaba.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_amb_file REFERENCE_GENOME_FASTA_AMB_FILE \
    --reference_genome_fasta_ann_file REFERENCE_GENOME_FASTA_ANN_FILE \
    --reference_genome_fasta_bwt_file REFERENCE_GENOME_FASTA_BWT_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --reference_genome_fasta_pac_file REFERENCE_GENOME_FASTA_PAC_FILE \
    --reference_genome_fasta_sa_file REFERENCE_GENOME_FASTA_SA_FILE \
    --params_svaba '"--hp --read-tracking"'
```

### Usage

```
workflow:
    1. Run Svaba.

usage: nexus run --nf-workflow paired-end_read_dna_variant_calling_svaba.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id',
                                            'tumor_bam_file',
                                            'tumor_bam_bai_file',
                                            'normal_bam_file',
                                            'normal_bam_bai_file'
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_amb_file   :   Reference genome FASTA.AMB file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.amb).
    --reference_genome_fasta_ann_file   :   Reference genome FASTA.ANN file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.ann).
    --reference_genome_fasta_bwt_file   :   Reference genome FASTA.BWT file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.bwt.2bit.64).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --reference_genome_fasta_pac_file   :   Reference genome FASTA.PAC file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.pac).
    --reference_genome_fasta_sa_file    :   Reference genome FASTA.PAC file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.sa).
    --params_svaba                      :   Svaba parameters (default: '"--hp --read-tracking"').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`

| Header              | Description                        |
|---------------------|------------------------------------|
| sample_id           | Sample ID                          |
| tumor_bam_file      | Full path to tumor `bam` file      |
| tumor_bam_bai_file  | Full path to tumor `bam.bai` file  |
| normal_bam_file     | Full path to normal `bam` file     |
| normal_bam_bai_file | Full path to normal `bam.bai` file |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--reference_genome_fasta_*_file`
* These index files must have been generated using `bwa` 

`--params_svaba`
* Refer to the [Svaba documentation](https://github.com/walaj/svaba).
* The following parameters for `svaba run` are already included in `nexus` module and should not be specified:
  * `--reference-genome`
  * `--id-string`
  * `--case-bam`
  * `--control-bam`
  * `--threads`
