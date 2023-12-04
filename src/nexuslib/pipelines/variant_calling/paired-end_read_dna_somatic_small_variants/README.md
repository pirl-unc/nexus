## paired_end_human_dna_somatic_small_variants.nf

Identifies small somatic DNA variants (SNVs and INDELs) with paired-end DNA sequencing 
(BAM) files using 
[GATK4-Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2), 
[Strelka2](https://github.com/Illumina/strelka), and
[DeepVariant](https://github.com/google/deepvariant).

### Inputs / Outputs

| Input(s)                  | Output(s)  |
|---------------------------|------------|
| `BAM` files | `VCF` files |

### Dependencies

* `gatk4-mutect2`
* `strelka2`
* `picard`
* `deepvariant`

### Usage

```
workflow:
    human:
        1.  Run gatk4-mutect2 (tumor and normal mode).
            Run gatk4 LearnReadOrientationModel.
            Run gatk4 GetPileupSummaries.
            Run gatk4 CalculateContamination.
            Run gatk4 FilterMutectCalls.
            Run picard MergeVcfs,
        2.  Run strelka2 (somatic mode).
        3.  Run deepvariant.
    non-human:
        1.  Run gatk4 mutect2 (tumor and normal mode).
            Run gatk4 FilterMutectCalls.
            Run picard MergeVcfs,
        2.  Run strelka2 (somatic mode).
        3.  Run deepvariant.

usage: nexus run --nf-workflow paired_end_dna_somatic_small_variants.nf [required] [optional] [--help]

required arguments:
    --samples_tsv_file                  :   TSV file with the following columns:
                                            'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file', 'tumor_sample_id', normal_sample_id'
    --output_dir                        :   Directory to which output files will be symlinked.

optional arguments:
    --reference_genome_fasta_file       :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --tools_list                        :   Tools to run (default: 'gatk4,strelka2,deepvariant').
    --is_human                          :   true if the samples are human. false otherwise (default: true).
    --python2                           :   python2 path (default: python2).
    --gatk4                             :   gatk4 path (default: gatk).
    --gatk4_mutect2_params              :   gatk4-mutect2 parameters (default:
                                            '"--germline-resource /datastore/lbcfs/collaborations/pirl/seqdata/references/af-only-gnomad.hg38.vcf
                                              --panel-of-normals /datastore/lbcfs/collaborations/pirl/seqdata/references/1000g_pon.hg38.vcf "').
                                            Note that the parameters need to be wrapped in quotes and a space at the end of the string is necessary.
    --gatk4_getpileupsummaries_params   :   gatk4 GetPileupSummaries parameters (default:
                                            '"-V /datastore/lbcfs/collaborations/pirl/seqdata/references/small_exac_common_3.hg38.vcf
                                              -L /datastore/lbcfs/collaborations/pirl/seqdata/references/small_exac_common_3.hg38.vcf "').
                                            Note that the parameters need to be wrapped in quotes and a space at the end of the string is necessary.
    --gatk4_chromosomes                 :   gatk4 chromosomes to parallelize (default:
                                            chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM).
    --picard                            :   Picard path (default:
                                            /datastore/lbcfs/collaborations/pirl/share/apps/picard/v3.1.1/picard.jar).
    --strelka2                          :   strelka2 configureStrelkaSomaticWorkflow.py path (default:
                                            /datastore/lbcfs/collaborations/pirl/share/apps/strelka2/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py).
    --strelka2_params                   :   strelka2 parameters (default: ' ').
                                            Note that the parameters need to be wrapped in quotes and a space at the end of the string is necessary.
    --containerization                  :   Containerization ('singularity' or 'docker'; default: 'singularity').
    --deepvariant_bin_path              :   DeepVariant bin path (e.g. '/opt/deepvariant/bin/run_deepvariant').
    --deepvariant_bin_version           :   DeepVariant bin version (default: '1.6.0').
    --deepvariant_input_path            :   DeepVariant input path (e.g. /datastore/).
    --deepvariant_output_path           :   DeepVariant output path (e.g. /datastore/).
    --deepvariant_model_type            :   DeepVariant --model_type parameter value (default: 'WGS').
    --delete_work_dir                   :   Delete work directory (default: false).
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