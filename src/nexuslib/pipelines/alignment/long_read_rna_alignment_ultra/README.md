## long_read_rna_alignment_ultra.nf

Aligns long RNA reads (Pacific Biosciences or Oxford Nanopore) to a reference genome using [uLTRA](https://github.com/ksahlin/ultra).

### Inputs / Outputs

| I/O    | Description                                                     |
|:-------|:----------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                                | 
| Output | MD-tagged and sorted `bam` and `bam.bai` files for each sample. |

### Dependencies

* `uLTRA`
* `samtools`

### Example

```
nexus run --nf-workflow long_read_rna_alignment_ultra.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --ultra_index ULTRA_INDEX \
    --params_ultra '"--isoseq"'
```

### Usage

```
workflow:
    1. Align reads (fastq.gz files) to a reference genome using uLTRA.
    2. Generate MD tags.
    3. Sort MD-tagged bam file.

usage: nexus run --nf-workflow long_read_rna_alignment_ultra.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'fastq_file'.
    --output_dir                        :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --reference_genome_fasta_fai_file   :   Reference genome FASTA.FAI file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa.fai).
    --ultra_index                       :   uLTRA index path (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/ultra/hg38_index/).
    --params_ultra                      :   uLTRA parameters (default: "--isoseq").
    --delete_work_dir                   :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow `config` file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--sample_tsv_file`

| Header     | Description                  |
| ---------- |------------------------------|
| sample_id  | Sample ID.                   |
| fastq_file | Full path to `fastq.gz` file |

`--reference_genome_fasta_file`
* Reference genome FASTA files can be found in /datastore/lbcfs/collaborations/pirl/seqdata/references/ on LBG.

`--ultra_index`
* An [uLTRA index](https://github.com/ksahlin/ultra) needs to be generated prior to running this workflow. 
* Prebuilt uLTRA indices are available in `/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/ultra/` on LBG.

`--params_ultra`
* Refer to the [uLTRA documentation](https://github.com/ksahlin/ultra).
* The following parameters for `minimap2` are already included in `nexus` module for `minimap2` and should not be specified:
  * `--t`
  * `--index`
  * `--prefix`