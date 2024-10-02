## paired-end_read_rna_alignment_star.nf

Aligns paired-end short DNA reads (Illumina) to a reference genome using [STAR](https://github.com/alexdobin/STAR).

### Inputs / Outputs

| I/O    | Description                                       |
|:-------|:--------------------------------------------------|
| Input  | `R1 and R2 fastq.gz` files for each sample.       | 
| Output | Sorted `bam` and `bam.bai` files for each sample. |

### Dependencies

* `STAR`
* `samtools`

### Example

```
nexus run --nf-workflow paired-end_read_rna_alignment_star.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --star_index STAR_INDEX \
    --params_star '"--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic"'
```

### Usage

```
workflow:
    1. Align paired-end reads to a reference genome index using STAR.

usage: nexus run --nf-workflow paired-end_read_rna_alignment_star.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --star_index                    :   Reference genome STAR index (default: /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/star/hg38_100bp_overhang/).
    --params_star                   :   STAR parameters (default: "--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic").
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                     |
| ------------ |---------------------------------|
| sample_id    | Sample ID.                      |
| fastq_file_1 | Full path to R1 `fastq.gz` file |
| fastq_file_2 | Full path to R2 `fastq.gz` file |

`--star_index`
* STAR index path. 
* Prebuilt indices are available in `/datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/star/` on LBG.

`--params_star`
* Refer to the [STAR documentation](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).
* The following parameters for `STAR` are already included in `nexus` module for `minimap2` and should not be specified:
  * `-t`
  * `-R`
