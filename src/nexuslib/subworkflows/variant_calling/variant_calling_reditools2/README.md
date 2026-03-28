## variant_calling_reditools2.nf

Identifies A-to-I RNA editing events in paired-end RNA and DNA BAM files using [reditools2](https://github.com/BioinfoUNIBA/REDItools2).

### Inputs / Outputs

| I/O    | Description                                            |
|:-------|:-------------------------------------------------------|
| Input  | RNA and DNA `bam` and `bam.bai` files for each sample. | 
| Output | `tsv` file for each sample.                            |

### Dependencies

* `reditools2`

### Example

```
nexus run --nf-workflow variant_calling_reditools2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genes_gtf_file REFERENCE_GENES_GTF_FILE \
    --params_reditools2 '""' \
    --params_reditools_annotatetable '"-s 4 -c 1,2,3 -n gencode"'
```

### Usage

```
workflow:
    1. Run REDItools2 reditools.py for RNA.
    2. Run REDItools2 reditools.py for DNA.
    3. Run annotate_with_DNA.py (identify RNA-specific A-to-I RNA editing events).
    4. Annotate variants.

usage: nexus run --nf-workflow variant_calling_reditools2.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id',
                                            'rna_bam_file',
                                            'rna_bam_bai_file',
                                            'dna_bam_file',
                                            'dna_bam_bai_file'
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --reference_genes_gtf_file          :   Reference genes GTF file.

optional arguments:
    --params_reditools2                 :   reditools.py parameters (default: '""').
                                            Note that the parameters need to be wrapped in quotes.
    --params_reditools_annotatetable    :   Reditools AnnotateTable.py parameters (default: '"-s 4 -c 1,2,3 -n gencode"').
                                            Note that the parameters need to be wrapped in quotes.
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                 |
|--------------|-----------------------------|
| sample_id    | Sample ID                   |
| bam_file     | Full path to `bam` file     |
| bam_bai_file | Full path to `bam.bai` file |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta` or `.fasta.gz`) file.

`--reference_genes_gtf_file`
* Reference genes GTF (`.gtf` or `.gtf.gz`) file.

`--params_reditools2`
* Refer to the [reditools documentation](https://github.com/BioinfoUNIBA/REDItools2).
* The following parameters for `reditools.py` are already included in `nexus` module for `reditools2` and should not be specified:
  * `-f`
  * `-o`
  * `-r`
  * `-S`
  * `-H`
  * `-B`
  * `--dna`

`--params_reditools_annotatetable`
* Refer to the [reditools documentation](https://github.com/BioinfoUNIBA/REDItools2).
* The following parameters for `AnnotateTable.py` (REDItools v1) are already included in `nexus` module for `reditools2` and should not be specified:
  * `-i`
  * `-a`
  * `-o`
