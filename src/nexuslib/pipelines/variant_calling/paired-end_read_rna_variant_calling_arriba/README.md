## paired-end_read_rna_variant_calling_arriba.nf

Identifies fusion genes in paired-read RNA fastq.gz files using [Arriba](https://arriba.readthedocs.io/en/latest/).

### Inputs / Outputs

| I/O    | Description                                             |
|:-------|:--------------------------------------------------------|
| Input  | RNA-seq `fastq.gz` files (R1 and R2) for each sample.   | 
| Output | `TSV` file for each sample.                             |

### Dependencies

* `Arriba`

### Example

```
nexus run --nf-workflow paired-end_read_rna_variant_calling_arriba.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --gtf_file GTF_FILE \
    --protein_domains_gff3_file PROTEIN_DOMAINS_GFF3_FILE \
    --params_arriba '"-S 3 -i chr*"'
```

### Usage

```
workflow:
    1. Run Arriba.

usage: nexus run --nf-workflow paired-end_read_rna_variant_calling_arriba.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --gtf_file                          :   Reference transcriptome GTF file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf).
    --protein_domains_gff3_file         :   Protein domains GFF3 file (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/arriba/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3).
    --params_arriba                     :   Arriba parameters (default: '"-S 3 -f blacklist -i chr*"').
                                            Note that the parameters need to be wrapped in quotes.
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                             |
|--------------|-----------------------------------------|
| sample_id    | Sample ID                               |
| fastq_file_1 | Full path to RNA-seq `fastq.gz` R1 file |
| fastq_file_2 | Full path to RNA-seq `fastq.gz` R2 file |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--reference_genome_fasta_fai_file`
* Reference genome FASTA.FAI files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--gtf_file`
* Reference transcriptome GTF files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--protein_domains_gff3_file`
* Reference protein domains GFF3 files for Arriba can be found in /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/arriba/arriba_v2.4.0/database on LBG.

`--params_arriba`
* Refer to the [Arriba documentation](https://arriba.readthedocs.io/en/latest/).
* The following parameters for `Arriba` are already included in `nexus` module for `delly call` and should not be specified:
  * `-x`
  * `-g`
  * `-a`
  * `-p`
  * `-o`