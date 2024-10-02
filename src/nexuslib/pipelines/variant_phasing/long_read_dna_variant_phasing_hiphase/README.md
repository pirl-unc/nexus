## long_read_dna_variant_phasing_hiphase.nf

Phases variants with long-read DNA BAM files using [HiPhase](https://github.com/PacificBiosciences/HiPhase).

### Inputs / Outputs

| I/O    | Description                                                                                                                    |
|:-------|:-------------------------------------------------------------------------------------------------------------------------------|
| Input  | `bam`, `bam.bai`, small variants `vcf.gz`, `vcf.gz.tbi`, and structural variants `vcf.gz`, `vcf.gz.tbi` files for each sample. | 
| Output | HiPhase output files (phased `vcf` and `bam`) for each sample.                                                                 |

### Dependencies

* `HiPhase`

### Example

```
nexus run --nf-workflow long_read_dna_variant_phasing_hiphase.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE
```

### Usage

```
workflow:
    1. Run HiPhase.

usage: nexus run --nf-workflow long_read_dna_variant_phasing_hiphase.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id',
                                            'bam_file',
                                            'bam_bai_file',
                                            'small_variants_vcf_gz_file',
                                            'small_variants_vcf_gz_tbi_file',
                                            'structural_variants_vcf_gz_file',
                                            'structural_variants_vcf_gz_tbi_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header                              | Description                                         |
|-------------------------------------|-----------------------------------------------------|
| sample_id                           | Sample ID.                                          |
| bam_file                            | Full path to `bam` file.                            |
| bam_bai_file                        | Full path to `bam.bai` file.                        |
| small_variants_vcf_gz_file          | Full path to small variants `vcf.gz` file.          |
| small_variants_vcf_gz_tbi_file      | Full path to small variants `vcf.gz.tbi` file.      |
| structural_variants_vcf_gz_file     | Full path to structural variants `vcf.gz` file.     |
| structural_variants_vcf_gz_tbi_file | Full path to structural variants `vcf.gz.tbi` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--reference_genome_fasta_fai_file`
* Reference genome FASTA.FAI files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.
