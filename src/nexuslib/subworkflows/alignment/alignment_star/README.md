## alignment_star.nf

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
nexus run --nf-workflow alignment_star.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --params_star_genomegenerate '"--genomeSAindexNbases 10 --sjdbOverhang 150"' \
    --params_star '"--genomeLoad NoSharedMemory --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic --chimSegmentMin 10 --chimOutType WithinBAM SoftClip Junctions --chimMultimapNmax 20 --chimOutJunctionFormat 0"'
```

### Usage

```
workflow:
    1. Index reference genome FASTA and GTF files.
    2. Align paired-end reads to a reference genome index using STAR.

usage: nexus run --nf-workflow alignment_star.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'fastq_file_1', 'fastq_file_2'.
    --output_dir                    :   Directory to which output files will be copied.
    --reference_genome_fasta_file   :   Reference genome FASTA file.
    --reference_genes_gtf_file      :   Reference genes GTF file.

optional arguments:
    --params_star_genomegenerate    :   STAR genomeGenerate parameter (default: "--genomeSAindexNbases 10 --sjdbOverhang 150").
                                        Note that the parameters need to be wrapped in quotes.
    --params_star                   :   STAR parameters (default: "--genomeLoad NoSharedMemory --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --twopassMode Basic --chimSegmentMin 10 --chimOutType WithinBAM SoftClip Junctions --chimMultimapNmax 20 --chimOutJunctionFormat 0").
                                        Note that the parameters need to be wrapped in quotes.
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

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta` or `.fasta.gz`) file.

`--reference_genes_gtf_file`
* * Reference genes GTF (`.gtf` or `.gtf.gz`) file.

`--params_star_genomegenerate`
* `STAR genomeGenerate` parameters (default: `--genomeSAindexNbases 10 --sjdbOverhang 150`). 
* Refer to the [STAR documentation](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).

`--params_star`
* Refer to the [STAR documentation](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).
* The following parameters for `STAR` are already included in `nexus` module for `minimap2` and should not be specified:
  * `-t`
  * `-R`
