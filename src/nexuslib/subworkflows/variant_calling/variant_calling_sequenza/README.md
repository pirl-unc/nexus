## variant_calling_sequenza.nf

Identifies somatic copy number alterations in paired-end DNA sequencing (BAM) files using [Sequenza](https://sequenza-utils.readthedocs.io/en/latest/).

### Inputs / Outputs

| I/O    | Description                                                |
|:-------|:-----------------------------------------------------------|
| Input  | Tumor and normal `bam` and `bam.bai` files for each sample. | 
| Output | Sequenza output files for each sample.                     |

### Dependencies

* `Sequenza`

### Example

```
nexus run --nf-workflow variant_calling_sequenza.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --assembly hg38 \
    --chromosomes 'chr1 chr2 chr3' \
    --params_sequenzautils_gcwiggle '"-w 50"' \
    --params_sequenza_bam2seqz '"-N 20 --qformat sanger"' \
    --params_sequenza_seqzbinning '"--window 50"'
```

### Usage

```
workflow:
    1.  Run sequenza-utils bam2seqz.
    2.  Merge seqz files into one seqz file.
    3.  Run sequenza-utils seqz_binning.
    4.  Run sequenza in R.

usage: nexus run --nf-workflow paired-end_read_dna_variant_calling_sequenza.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id',
                                            'tumor_bam_file',
                                            'tumor_bam_bai_file',
                                            'normal_bam_file',
                                            'normal_bam_bai_file'
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --assembly                          :   Assembly ('hg19' or 'hg38').

optional arguments:
    --chromosomes                       :   List of chromosomes (default: '"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 ch22 chrX chrY"').
                                            Note that the parameters need to be wrapped in quotes.
    --params_sequenzautils_gcwiggle     :   sequenza-utils gc_wiggle parameters (default: '"-w 50"').
                                            Note that the parameters need to be wrapped in quotes.
    --params_sequenza_bam2seqz          :   Sequenza 'bam2seqz' parameters (default: '"-N 20 --qformat sanger"').
                                            Note that the parameters need to be wrapped in quotes.
    --params_sequenza_seqzbinning       :   Sequenza 'seqz_binning' parameters (default: '"--window 50"').
                                            Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header              | Description                           |
|---------------------|---------------------------------------|
| sample_id           | Sample ID.                            |
| tumor_bam_file      | Full path to tumor `bam` file.        |
| tumor_bam_bai_file  | Full path to tumor `bam.bai` file.    |
| normal_bam_file     | Full path to normal `bam` file.       |
| normal_bam_bai_file | Full path to normal `bam.bai` file.   |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--assembly`
* Genome assembly.

* `--params_sequenzautils_gcwiggle`
* Refer to the [Sequenza documentation](https://sequenza-utils.readthedocs.io/en/latest/).
* The following parameters for `Sequenza` are already included in `nexus` module for `sequenza-utils gc_wiggle` and should not be specified:
  * `-f`
  * `-o`

`--params_sequenza_bam2seqz`
* Refer to the [Sequenza documentation](https://sequenza-utils.readthedocs.io/en/latest/).
* The following parameters for `Sequenza` are already included in `nexus` module for `sequenza-utils bam2seqz` and should not be specified:
  * `-t`
  * `-n`
  * `--fasta`
  * `-gc`
  * `--chromosome`
  * `--parallel`
  * `-o`

`--params_sequenza_seqzbinning`
* Refer to the [Sequenza documentation](https://sequenza-utils.readthedocs.io/en/latest/).
* The following parameters for `Sequenza` are already included in `nexus` module for `sequenza-utils seqz_binning` and should not be specified:
  * `--seqz`
  * `-o`
