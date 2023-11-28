# Nexus

**NEX**tflow's **U**ltimate **S**treamliner

Command-line interface (CLI) for conveniently running bioinformatics 
workflows written in nextflow.

[![build](https://github.com/pirl-unc/nexus/actions/workflows/main.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/main.yml)

## 01. Installation

```
pip install . --verbose
```

## 02. Dependencies

```
conda install pandas
conda install nextflow
conda install samtools==1.18
conda install minimap2==2.22
conda install sniffles==2.2
conda install svim==2.0.0
conda install pbsv==2.9.0
conda install cutesv==2.1.0
pip install ultra-bioinformatics
```

## 03. Usage

To view all available workflows
```
nexus avail
```

To run a workflow
```
nexus run [-h] [--nf-workflow NF_NEXTFLOW] [workflow-specific parameters]
```

To view workflow-specific parameters
```
nexus run --nf-workflow <workflow> --help
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
    -c nextflow.config \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file FASTA_FILE \
    --minimap2 minimap2 \
    --minimap2_params "'-ax map-hifi --cs --eqx -Y -L'" \
    --samtools samtools \
    --platform_tag pacbio \
```

`nextflow.config` can be downloaded from [here](/configs/) for 
local computers (laptop) and HPC (slurm).

`SAMPLES_TSV_FILE` is a tab-delimited file:

| Header     | Description                  |
| ---------- |------------------------------|
| sample_id  | Sample ID                    |
| fastq_file | Full path to `fastq.gz` file |

For more on this particular workflow, check out [here](/src/nexuslib/pipelines/alignment/long_read_alignment_minimap2/).

## 04. Available workflows

Here is a list of all available workflows for `v0.0.1`

| Category       | Workflow                                                                                                                                    |
|:---------------|:--------------------------------------------------------------------------------------------------------------------------------------------|
| Alignment      | [long_read_alignment_minimap2.nf](/src/nexuslib/pipelines/alignment/long_read_alignment_minimap2/)                                          |
| Alignment      | [long_read_rna_alignment_ultra.nf](/src/nexuslib/pipelines/alignment/long_read_rna_alignment_ultra)                                         |
| Alignment      | [paired-end_read_dna_alignment_bwa-mem2.nf](/src/nexuslib/pipelines/alignment/paired-end_read_dna_alignment_bwa-mem2/)                      |
| Quantification | [paired-end_read_rna_quantification_salmon_mapping_mode.nf](/src/nexuslib/pipelines/quantification/paired-end_read_rna_quantification_salmon-mapping/) |
