# Nexus

**NEX**tflow's **U**ltimate **S**treamliner

Command-line interface (CLI) for conveniently running bioinformatics 
workflows written in nextflow.

[![build](https://github.com/pirl-unc/nexus/actions/workflows/main.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/main.yml)

## 01. Installation

Download the [latest stable release](https://github.com/pirl-unc/nexus/releases/download/v0.0.1/nexus-0.0.1.tar.gz)

```
conda create -n nexus python=3.10
conda activate nexus
pip install nexus-<version>.tar.gz --verbose
```

## 02. Dependencies

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install nextflow
conda install pandas
conda install samtools==1.18
conda install minimap2==2.22
conda install bwa-mem2==2.2.1
conda install gatk4==4.4.0.0
conda install sniffles==2.2
conda install svim==2.0.0
conda install pbsv==2.9.0
conda install cutesv==2.1.0
conda install delly==1.1.8
conda install bcftools==1.18
conda install sambamba==1.0
conda install samblaster==0.1.26
pip install ultra-bioinformatics

conda create -n py27 python=2.7
conda activate py27
conda install lumpy-sv==0.3.1
conda deactivate
```

## 03. Usage

To view all available workflows
```
nexus avail
```

To run a workflow
```
nexus run [-h] [--nf-workflow NF_WORKFLOW] [workflow-specific parameters]
```

To view workflow-specific parameters
```
nexus run --nf-workflow [NF_WORKFLOW] --help
```

### Example

Long-read whole-genome sequencing `fastq.gz` alignment

| I/O    | Description                                                                  |
|:-------|:-----------------------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.<br/>                                        | 
| Ouptut | MD-tagged and sorted `bam` and `bam.bai` files for each sample. |

```
nexus run \
    --nf-workflow long_read_alignment_minimap2.nf \
    -c nextflow_slurm.config \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR
```

`nextflow_slurm.config` can be downloaded from [here](/configs/nextflow_slurm.config).

`SAMPLES_TSV_FILE` is a tab-delimited file:

| Header     | Description                  |
| ---------- |------------------------------|
| sample_id  | Sample ID                    |
| fastq_file | Full path to `fastq.gz` file |

For more on this particular workflow, check out [here](/src/nexuslib/pipelines/alignment/long_read_alignment_minimap2/).

## 04. Available workflows

Here is a list of all available workflows for `v0.0.2`

| Category       | Workflow                                                                                                               |
|:---------------|:-----------------------------------------------------------------------------------------------------------------------|
| Alignment      | [long_read_alignment_minimap2.nf](/src/nexuslib/pipelines/alignment/long_read_alignment_minimap2/)                     |
| Alignment      | [long_read_rna_alignment_ultra.nf](/src/nexuslib/pipelines/alignment/long_read_rna_alignment_ultra)                    |
| Alignment      | [paired-end_read_dna_alignment_bwa-mem2.nf](/src/nexuslib/pipelines/alignment/paired-end_read_dna_alignment_bwa-mem2/) |
| Utilities       | [sequencing_coverage.nf](/src/nexuslib/pipelines/utilities/sequencing_coverage/) |
| Variant calling | [long_read_dna_small_variants.nf](/src/nexuslib/pipelines/variant_calling/long_read_dna_small_variants/)               |
| Variant calling | [long_read_dna_structural_variants.nf](/src/nexuslib/pipelines/variant_calling/long_read_dna_structural_variants/)     |
| Variant calling | [paired-end_read_dna_somatic_small_variants.nf](/src/nexuslib/pipelines/variant_calling/paired-end_read_dna_somatic_small_variants/) |
| Variant calling | [paired-end_read_dna_somatic_structural_variants.nf](/src/nexuslib/pipelines/variant_calling/paired-end_read_dna_somatic_structural_variants/) |