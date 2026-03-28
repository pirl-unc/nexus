# Nexus

**NEX**tflow's **U**ltimate **S**treamliner

A command-line interface (CLI) for running bioinformatics workflows written in Nextflow.

[![CI](https://github.com/pirl-unc/nexus/actions/workflows/ci.yml/badge.svg)](https://github.com/pirl-unc/nexus/actions/workflows/ci.yml)

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

## 04. License

Licensed under the Apache License, Version 2.0.
