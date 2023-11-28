## paired_end_human_dna_somatic_small_variants.nf

Identifies small DNA variants (SNVs and INDELs) with paired-end DNA sequencing 
(BAM) files using 
[GATK4-Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2), 
[Strelka2](https://github.com/Illumina/strelka).

### Inputs / Outputs

| Input(s)                  | Output(s)  |
|---------------------------|------------|
| `BAM` files | `VCF` files |

### Dependencies

* `gatk4-mutect2`
* `strelka2`

### Usage

```shell
nextflow run paired_end_human_dna_somatic_small_variants.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file', 'tumor_sample_id', normal_sample_id'
    --output_dir                        :   Directory to which output files will be symlinked.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --gatk4                             :   GATK4 path.
    --gatk4_mutect2_params              :   GATK4-Mutect2 parameters (e.g. "--germline-resource /<path>/af-only-gnomad.hg38.vcf --panel-of-normals /<path>/1000g_pon.hg38.vcf.gz").
    --gatk4_getpileupsummaries_params   :   GATK4-GetPileupSummaries parameters (e.g. "-V /<path>/small_exac_common_3.hg38.vcf -L /<path>/small_exac_common_3.hg38.vcf").
    --picard                            :   Picard path.
    --strelka2                          :   Strelka2 path (bin path).
    --strelka2_params                   :   Strelka2 parameters.
    --chromosomes                       :   Chromosomes to parallelize using GATK4 (separated by comma; e.g. 'chr1,chr2,chr3').
```

### Parameters

`--sample_tsv_file`

| Header       | Description                 |
|--------------|-----------------------------|
| sample_id    | Sample ID                   |
| bam_file     | Full path to `BAM` file     |
| bam_bai_file | Full path to `BAM.BAI` file |
