## paired_end_human_dna_somatic_small_variants.nf

Identifies structural DNA variants with paired-end DNA sequencing 
(BAM) files using 
[Delly2](https://github.com/dellytools/delly), 
[Lumpy](https://github.com/arq5x/lumpy-sv).

### Inputs / Outputs

| Input(s)                  | Output(s)  |
|---------------------------|------------|
| `BAM` files | `VCF` files |

### Dependencies

* `delly`
* `lumpyexpress`

### Usage

```
workflow:
    1. Run delly2 (tumor and normal mode).
    2. Run lumpyexpress (tumor and normal mode).

usage: nexus run --nf-workflow paired_end_dna_somatic_structural_variants.nf [required] [optional] [--help]

required arguments:
    -c                                          :   Nextflow .config file.
    -w                                          :   Nextflow work directory path.
    --samples_tsv_file                          :   TSV file with the following columns:
                                                    'sample_id',
                                                    'tumor_bam_file',
                                                    'tumor_bam_bai_file',
                                                    'normal_bam_file',
                                                    'normal_bam_bai_file',
                                                    'tumor_sample_id',
                                                    'normal_sample_id'
    --output_dir                                :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file               :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --tools_list                                :   Tools to run (default: 'delly2,lumpyexpress').
    --delly2                                    :   delly2 path (default: delly).
    --delly2_params                             :   delly2 'call' parameters (default:
                                                    '"--exclude /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/delly2/human.hg38.excl.tsv --map-qual 20 "').
                                                    Note that the parameters need to be wrapped in quotes
                                                    and a space at the end of the string is necessary.
    --bcftools                                  :   bcftools path (default: bcftools).
    --python2                                   :   python2 path (default: python2).
    --lumpyexpress                              :   lumpyexpress path (default: lumpyexpress).
    --lumpyexpress_config_file                  :   lumpyexpress config file (default: lumpyexpress.config).
    --lumpy_extract_split_reads_script_file     :   lumpy 'extractSplitReads_BwaMem' file (default:
                                                    /datastore/lbcfs/collaborations/pirl/share/apps/lumpy/lumpy-sv-0.3.1/scripts/extractSplitReads_BwaMem).
    --samtools                                  :   samtools path (default: samtools).
    --delete_work_dir                           :   Delete work directory (default: false).
```

### Parameters

`--sample_tsv_file`

| Header              | Description                           |
|---------------------|---------------------------------------|
| sample_id           | Sample ID.                            |
| tumor_bam_file      | Full path to tumor `bam` file.        |
| tumor_bam_bai_file  | Full path to tumor `bam.bai` file.    |
| normal_bam_file     | Full path to normal `bam` file.       |
| normal_bam_bai_file | Full path to normal `bam.bai` file.   |
| tumor_sample_id     | Tumor sample ID.                      |
| normal_sample_id    | Normal sample ID.                     |
