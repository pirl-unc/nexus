# Nexus

**NEX**tflow's **U**ltimate **S**treamliner

A command-line interface (CLI) for running bioinformatics workflows written in Nextflow.

## 01. Installation

Download the latest stable release [here](https://github.com/pirl-unc/nexus/releases).

```bash
conda create -n nexus python=3.10
conda activate nexus
conda install bioconda::nextflow==24.10.0
pip install nexus-<version>.tar.gz --verbose
```

## 02. Dependencies

* Python (>=3.10)
* Nextflow (>=24.10.0)

## 03. Usage

### View all available workflows
```bash
nexus avail
```

### Run a subworkflow or workflow
```bash
nexus run --nf-workflow <WORKFLOW.nf> [workflow-specific parameters]
```

### View workflow-specific parameters
```bash
nexus run --nf-workflow <WORKFLOW.nf> --help
```

## 04. CI Status

### Subworkflows

| Category | Status |
|:---------|:-------|
| Alignment | [![alignment](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_alignment.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_alignment.yml) |
| Antigen Prediction | [![antigen_prediction](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_antigen_prediction.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_antigen_prediction.yml) |
| Assembly | [![assembly](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_assembly.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_assembly.yml) |
| Haplotyping | [![haplotyping](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_haplotyping.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_haplotyping.yml) |
| MHC Typing | [![mhc_typing](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_hla_typing.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_hla_typing.yml) |
| Isoform Characterization | [![isoform_characterization](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_isoform_characterization.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_isoform_characterization.yml) |
| Peptide Prediction | [![peptide_prediction](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_peptide_prediction.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_peptide_prediction.yml) |
| Quantification | [![quantification](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_quantification.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_quantification.yml) |
| Read Error Correction | [![read_error_correction](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_read_error_correction.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_read_error_correction.yml) |
| Sequencing Simulation | [![sequencing_simulation](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_sequencing_simulation.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_sequencing_simulation.yml) |
| Utilities | [![utilities](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_utilities.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_utilities.yml) |
| Variant Annotation | [![variant_annotation](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_variant_annotation.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_variant_annotation.yml) |
| Variant Calling | [![variant_calling](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_variant_calling.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_variant_calling.yml) |
| Variant Phasing | [![variant_phasing](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_variant_phasing.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/subworkflow_variant_phasing.yml) |

### Workflows

| Category | Status |
|:---------|:-------|
| MHC Typing | [![mhc_typing](https://github.com/pirl-unc/nexus/actions/workflows/workflow_hla_typing.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/workflow_hla_typing.yml) |
| Isoform Characterization | [![isoform_characterization](https://github.com/pirl-unc/nexus/actions/workflows/workflow_isoform_characterization.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/workflow_isoform_characterization.yml) |
| Quantification | [![quantification](https://github.com/pirl-unc/nexus/actions/workflows/workflow_quantification.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/workflow_quantification.yml) |
| Variant Calling | [![variant_calling](https://github.com/pirl-unc/nexus/actions/workflows/workflow_variant_calling.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/workflow_variant_calling.yml) |


## 05. License

Licensed under the Apache License, Version 2.0.
