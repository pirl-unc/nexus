## variant_calling_arriba.nf

Identifies fusion genes in paired-read RNA fastq.gz files using [Arriba](https://arriba.readthedocs.io/en/latest/).

### Inputs / Outputs

| I/O    | Description                                       |
|:-------|:--------------------------------------------------|
| Input  | RNA-seq `BAM` and `BAM.BAI` files for each sample. | 
| Output | `TSV` file for each sample.                       |

### Dependencies

* `Arriba`

### Example

```
nexus run --nf-workflow variant_calling_arriba.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --protein_domains_gff3_file PROTEIN_DOMAINS_GFF3_FILE \
    --params_arriba '"-S 3 -f blacklist -i chr*"'
```

### Usage

```
workflow:
    1. Align paired-end reads to a reference genome index using STAR.
    2. Run Arriba.

usage: nexus run --nf-workflow variant_calling_arriba.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --reference_genes_gtf_file          :   Reference genes GTF file.
    --protein_domains_gff3_file         :   Protein domains GFF3 file.

optional arguments:
    --params_arriba                     :   Arriba parameters (default: '"-S 3 -f blacklist -i chr*"').
                                            Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                           |
|--------------|---------------------------------------|
| sample_id    | Sample ID                             |
| bam_file     | Full path to RNA-seq `BAM` file       |
| bam_bai_file | Full path to RNA-seq `BAM.BAI` file   |

Run STAR with following parameters:

```
--genomeLoad NoSharedMemory
--readFilesCommand zcat 
--outSAMtype BAM SortedByCoordinate 
--outSAMunmapped Within 
--outSAMattributes Standard 
--twopassMode Basic 
--chimSegmentMin 10 
--chimOutType WithinBAM SoftClip Junctions
--chimMultimapNmax 20 
--chimOutJunctionFormat 1
```

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta` or `.fasta.gz`) file.

`--reference_genes_gtf_file`
* Reference genes GTF (`.gtf` or `.gtf.gz`) file.

`--protein_domains_gff3_file`
* Reference protein domains GFF3 file for Arriba. This file can be downloaded from https://github.com/suhrig/arriba/releases.

`--params_arriba`
* Refer to the [Arriba documentation](https://github.com/suhrig/arriba/wiki/01-Home).
* The following parameters for `Arriba` are already included in `nexus` module for `arriba` and should not be specified:
  * `-x`
  * `-g`
  * `-a`
  * `-p`
  * `-o`
