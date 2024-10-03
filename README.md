# Nexus

**NEX**tflow's **U**ltimate **S**treamliner

Command-line interface (CLI) for running bioinformatics workflows written in nextflow.

* [![alignment](https://github.com/pirl-unc/nexus/actions/workflows/alignment.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/alignment.yml)
* [![antigen_presentation_prediction](https://github.com/pirl-unc/nexus/actions/workflows/antigen_presentation_prediction.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/antigen_presentation_prediction.yml)
* [![assembly](https://github.com/pirl-unc/nexus/actions/workflows/assembly.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/assembly.yml)
* [![haplotyping](https://github.com/pirl-unc/nexus/actions/workflows/haplotyping.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/haplotyping.yml)
* [![hla_typing](https://github.com/pirl-unc/nexus/actions/workflows/hla_typing.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/hla_typing.yml)
* [![novel_isoform_discovery](https://github.com/pirl-unc/nexus/actions/workflows/novel_isoform_discovery.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/novel_isoform_discovery.yml)
* [![quantification](https://github.com/pirl-unc/nexus/actions/workflows/quantification.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/quantification.yml)
* [![read_error_correction](https://github.com/pirl-unc/nexus/actions/workflows/read_error_correction.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/read_error_correction.yml)
* [![utilities](https://github.com/pirl-unc/nexus/actions/workflows/utilities.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/utilities.yml)
* [![variant_calling_dna_long_reads](https://github.com/pirl-unc/nexus/actions/workflows/variant_calling_dna_long_reads.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/variant_calling_dna_long_reads.yml)
* [![variant_calling_dna_paired_end_reads](https://github.com/pirl-unc/nexus/actions/workflows/variant_calling_dna_paired_end_reads.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/variant_calling_dna_paired_end_reads.yml)
* [![variant_calling_rna_long_reads](https://github.com/pirl-unc/nexus/actions/workflows/variant_calling_rna_long_reads.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/variant_calling_rna_long_reads.yml)
* [![variant_calling_rna_paired_end_reads](https://github.com/pirl-unc/nexus/actions/workflows/variant_calling_rna_paired_end_reads.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/variant_calling_rna_paired_end_reads.yml)
* [![variant_phasing_dna](https://github.com/pirl-unc/nexus/actions/workflows/variant_phasing_dna.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/variant_phasing_dna.yml)

## 01. Installation

Download the latest stable release [here](https://github.com/pirl-unc/nexus/releases)

```
conda create -n nexus python=3.10
conda activate nexus
conda install bioconda::nextflow==23.10.0
pip install nexus-<version>.tar.gz --verbose
```

## 02. Dependencies

* python (>=3.10)
* nextflow (DSL 2)

## 03. Usage

### To view all available workflows
```
nexus avail
```

### To run a workflow
```
nexus run [-h] [--nf-workflow NF_WORKFLOW] [workflow-specific parameters]
```

### To view workflow-specific parameters
```
nexus run --nf-workflow [NF_WORKFLOW] --help
```

## 04. Example

Long-read whole-genome sequencing `fastq.gz` alignment to a reference genome.

| I/O    | Description                                                                  |
|:-------|:-----------------------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.<br/>                                        | 
| Ouptut | MD-tagged and sorted `bam` and `bam.bai` files for each sample. |

```
nexus run --nf-workflow long_read_alignment_minimap2.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --reference_genome_fasta_fai_file REFERENCE_GENOME_FASTA_FAI_FILE \
    --params_minimap2 '"-ax map-hifi --cs --eqx -Y -L"'
```

`NEXTFLOW_CONFIG_FILE` can be downloaded from [here](/nextflow/).

`SAMPLES_TSV_FILE` is a tab-delimited file:

| Header     | Description                  |
| ---------- |------------------------------|
| sample_id  | Sample ID                    |
| fastq_file | Full path to `fastq.gz` file |

For more on this particular workflow, check out [here](/src/nexuslib/pipelines/alignment/long_read_alignment_minimap2/).

## 05. Documentation for Available Workflows

A list of links to documentation for all available workflows in the latest version is provided below:

| Category                        | Workflow                                                                                                                                                 |
|:--------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------|
| Alignment                       | [**long_read_alignment_minimap2.nf**](/src/nexuslib/pipelines/alignment/long_read_alignment_minimap2/)                                                   |
| Alignment                       | [**long_read_rna_alignment_ultra.nf**](/src/nexuslib/pipelines/alignment/long_read_rna_alignment_ultra)                                                  |
| Alignment                       | [**paired-end_read_dna_alignment_bwa-mem2.nf**](/src/nexuslib/pipelines/alignment/paired-end_read_dna_alignment_bwa-mem2/)                               |
| Alignment                       | [**paired-end_read_rna_alignment_star.nf**](/src/nexuslib/pipelines/alignment/paired-end_read_rna_alignment_star/)                                       |
| Antigen presentation prediction | [**antigen_presentation_prediction_mhcflurry2.nf**](/src/nexuslib/pipelines/antigen_presentation_prediction/antigen_presentation_prediction_mhcflurry2/) |
| Assembly                        | [**transcriptome_assembly_rnabloom2.nf**](/src/nexuslib/pipelines/assembly/transcriptome_assembly_rnabloom2/)                                            | 
| Assembly                        | [**transcriptome_assembly_stringtie2.nf**](/src/nexuslib/pipelines/assembly/transcriptome_assembly_stringtie2/)                                          |
| Haplotyping                     | [**long_read_dna_haplotyping_whatshap.nf**](/src/nexuslib/pipelines/haplotyping/long_read_dna_haplotyping_whatshap/)                                     | 
| HLA typing                      | [**paired-end_read_rna_hla_typing_arcashla.nf**](/src/nexuslib/pipelines/hla_typing/paired-end_read_rna_hla_typing_arcashla/)                            |
| HLA typing                      | [**paired-end_read_rna_hla_typing_hlaprofiler.nf**](/src/nexuslib/pipelines/hla_typing/paired-end_read_rna_hla_typing_hlaprofiler/)                      |
| HLA typing                      | [**paired-end_read_rna_hla_typing_seq2hla.nf**](/src/nexuslib/pipelines/hla_typing/paired-end_read_rna_hla_typing_seq2hla/)                              |
| Novel isoform discovery         | [**novel_isoform_discovery_flair.nf**](/src/nexuslib/pipelines/novel_isoform_discovery/novel_isoform_discovery_flair/)                                   |
| Novel isoform discovery         | [**novel_isoform_discovery_isoquant.nf**](/src/nexuslib/pipelines/novel_isoform_discovery/novel_isoform_discovery_isoquant/)                             |
| Novel isoform discovery         | [**novel_isoform_discovery_isoseq.nf**](/src/nexuslib/pipelines/novel_isoform_discovery/novel_isoform_discovery_isoseq/)                                 |
| Quantification (RNA)            | [**long_read_rna_quantification_bambu.nf**](/src/nexuslib/pipelines/quantification/long_read_rna_quantification_bambu/)                                  |
| Quantification (RNA)            | [**long_read_rna_quantification_kallisto.nf**](/src/nexuslib/pipelines/quantification/long_read_rna_quantification_kallisto/)                            |
| Quantification (RNA)            | [**long_read_rna_quantification_liqa.nf**](/src/nexuslib/pipelines/quantification/long_read_rna_quantification_liqa/)                                    |
| Quantification (RNA)            | [**long_read_rna_quantification_oarfish.nf**](/src/nexuslib/pipelines/quantification/long_read_rna_quantification_oarfish/)                              |
| Quantification (RNA)            | [**paired-end_read_rna_quantification_salmon_mapping.nf**](/src/nexuslib/pipelines/quantification/paired-end_read_rna_quantification_salmon_mapping/)    |
| Quantification (RNA)            | [**paired-end_read_rna_quantification_kallisto.nf**](/src/nexuslib/pipelines/quantification/paired-end_read_rna_quantification_kallisto/)                |
| Read error correction           | [**long_read_error_correction_ratatosk.nf**](/src/nexuslib/pipelines/read_error_correction/long_read_error_correction_ratatosk/)                         |
| Utilities                       | [**fastq_to_unaligned_bam.nf**](/src/nexuslib/pipelines/utilities/fastq_to_unaligned_bam/)                                                               |
| Utilities                       | [**fastqc.nf**](/src/nexuslib/pipelines/utilities/fastqc/)                                                                                               |
| Utilities                       | [**sequencing_coverage.nf**](/src/nexuslib/pipelines/utilities/sequencing_coverage/)                                                                     |
| Variant calling (DNA)           | [**long_read_dna_variant_calling_clairs.nf**](/src/nexuslib/pipelines/variant_calling/long_read_dna_variant_calling_clairs/)                             |
| Variant calling (DNA)           | [**long_read_dna_variant_calling_cutesv.nf**](/src/nexuslib/pipelines/variant_calling/long_read_dna_variant_calling_cutesv/)                             |
| Variant calling (DNA)           | [**long_read_dna_variant_calling_deepvariant.nf**](/src/nexuslib/pipelines/variant_calling/long_read_dna_variant_calling_deepvariant/)                   |
| Variant calling (DNA)           | [**long_read_dna_variant_calling_dysgu.nf**](/src/nexuslib/pipelines/variant_calling/long_read_dna_variant_calling_dysgu/)                               |
| Variant calling (DNA)           | [**long_read_dna_variant_calling_nanomonsv.nf**](/src/nexuslib/pipelines/variant_calling/long_read_dna_variant_calling_nanomonsv/)                       |
| Variant calling (DNA)           | [**long_read_dna_variant_calling_pbsv.nf**](/src/nexuslib/pipelines/variant_calling/long_read_dna_variant_calling_pbsv/)                                 |
| Variant calling (DNA)           | [**long_read_dna_variant_calling_savana.nf**](/src/nexuslib/pipelines/variant_calling/long_read_dna_variant_calling_savana/)                             |
| Variant calling (DNA)           | [**long_read_dna_variant_calling_severus.nf**](/src/nexuslib/pipelines/variant_calling/long_read_dna_variant_calling_severus/)                           |
| Variant calling (DNA)           | [**long_read_dna_variant_calling_sniffles2.nf**](/src/nexuslib/pipelines/variant_calling/long_read_dna_variant_calling_sniffles2/)                       |
| Variant calling (DNA)           | [**long_read_dna_variant_calling_svim.nf**](/src/nexuslib/pipelines/variant_calling/long_read_dna_variant_calling_svim/)                                 |
| Variant calling (DNA)           | [**long_read_dna_variant_calling_svisionpro.nf**](/src/nexuslib/pipelines/variant_calling/long_read_dna_variant_calling_svisionpro/)                     |
| Variant calling (DNA)           | [**paired-end_read_dna_variant_calling_delly2.nf**](/src/nexuslib/pipelines/variant_calling/paired-end_read_dna_variant_calling_delly2/)                 |
| Variant calling (DNA)           | [**paired-end_read_dna_variant_calling_dysgu.nf**](/src/nexuslib/pipelines/variant_calling/paired-end_read_dna_variant_calling_dysgu/)                   |
| Variant calling (DNA)           | [**paired-end_read_dna_variant_calling_gridss.nf**](/src/nexuslib/pipelines/variant_calling/paired-end_read_dna_variant_calling_gridss/)                 |
| Variant calling (DNA)           | [**paired-end_read_dna_variant_calling_lumpy.nf**](/src/nexuslib/pipelines/variant_calling/paired-end_read_dna_variant_calling_lumpy/)                   |
| Variant calling (DNA)           | [**paired-end_read_dna_variant_calling_manta.nf**](/src/nexuslib/pipelines/variant_calling/paired-end_read_dna_variant_calling_manta/)                   |
| Variant calling (DNA)           | [**paired-end_read_dna_variant_calling_mutect2.nf**](/src/nexuslib/pipelines/variant_calling/paired-end_read_dna_variant_calling_mutect2/)               |
| Variant calling (DNA)           | [**paired-end_read_dna_variant_calling_sequenza.nf**](/src/nexuslib/pipelines/variant_calling/paired-end_read_dna_variant_calling_sequenza/)             |
| Variant calling (DNA)           | [**paired-end_read_dna_variant_calling_strelka2.nf**](/src/nexuslib/pipelines/variant_calling/paired-end_read_dna_variant_calling_strelka2/)             |
| Variant calling (DNA)           | [**paired-end_read_dna_variant_calling_svaba.nf**](/src/nexuslib/pipelines/variant_calling/paired-end_read_dna_variant_calling_svaba/)                   |
| Variant calling (RNA)           | [**long_read_rna_variant_calling_de_souza.nf**](/src/nexuslib/pipelines/variant_calling/long_read_rna_variant_calling_de_souza/)                         |
| Variant calling (RNA)           | [**paired-end_read_rna_variant_calling_arriba.nf**](/src/nexuslib/pipelines/variant_calling/paired-end_read_rna_variant_calling_arriba/)                 |
| Variant phasing (DNA)           | [**long_read_dna_variant_phasing_hiphase.nf**](/src/nexuslib/pipelines/variant_phasing/long_read_dna_variant_phasing_hiphase/)                           |
| Variant phasing (DNA)           | [**long_read_dna_variant_phasing_whatshap.nf**](/src/nexuslib/pipelines/variant_phasing/long_read_dna_variant_phasing_whatshap/)                         |

