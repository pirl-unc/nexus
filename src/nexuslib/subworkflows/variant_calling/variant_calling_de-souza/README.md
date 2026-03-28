## variant_calling_de-souza.nf

Identifies small RNA variants (SNVs and INDELs) in long-read RNA BAM files using [lrRNAseqVariantCalling](https://github.com/vladimirsouza/lrRNAseqVariantCalling).

### Inputs / Outputs

| I/O    | Description                                                       |
|:-------|:------------------------------------------------------------------|
| Input  | `fastq.gz` file for each sample.                                  | 
| Output | Flag corrected `bam`, `bam.bai`, and `vcf` files for each sample. |

### Dependencies

* `lrRNAseqVariantCalling (de Souza et al., Genome Biology 2023)`

### Example

```
nexus run --nf-workflow variant_calling_de-souza.nf \
    -c NEXTFLOW_CONFIG_FILE \
    -w WORK_DIR \
    --samples_tsv_file SAMPLES_TSV_FILE \
    --output_dir OUTPUT_DIR \
    --reference_genome_fasta_file REFERENCE_GENOME_FASTA_FILE \
    --deepvariant_containerization {singularity,docker} \
    --deepvariant_model_type PACBIO \
    --deepvariant_bin_path /opt/deepvariant/bin/run_deepvariant\
    --deepvariant_bin_version 1.6.0 \
    --deepvariant_input_path INPUT_PATH \
    --deepvariant_output_path OUTPUT_PATH
```

### Usage

```
workflow:
    1. Align reads to the reference using minimap2.
    2. Filter, sort, and index using samtools.
    3. Run GATK4 SplitNCigarReads.
    4. Perform flag correction.
    5. Index using samtools.
    6. Run DeepVariant.

usage: nexus run --nf-workflow variant_calling_de-souza.nf [required] [optional] [--help]

required arguments:
    -c                                  :   Nextflow .config file.
    -w                                  :   Nextflow work directory path.
    --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'fastq_file'.
    --output_dir                        :   Directory to which output files will be copied.
    --reference_genome_fasta_file       :   Reference genome FASTA file.
    --deepvariant_input_path            :   DeepVariant input path.
    --deepvariant_output_path           :   DeepVariant output path.

optional arguments:
    --params_minimap2                   :   Minimap2 parameters (default: '"-ax splice -uf -C5 --secondary=no"').
                                            Note that the parameters need to be wrapped in quotes.
    --params_samtools_view              :   Samtools view parameters (default: '"-F 2308"').
                                            Note that the parameters need to be wrapped in quotes.
    --platform_tag                      :   Platform tag (default: 'unknown').
    --platform_unit_tag                 :   Platform unit tag (default: 'unknown').
    --library_tag                       :   Library tag (default: 'unknown').
    --deepvariant_containerization      :   DeepVariant containerization ('singularity' or 'docker'; default: 'singularity').
    --deepvariant_model_type            :   DeepVariant --model_type parameter value (default: 'PACBIO').
    --deepvariant_bin_path              :   DeepVariant bin path (default: '/opt/deepvariant/bin/run_deepvariant').
    --deepvariant_bin_version           :   DeepVariant bin version (default: '1.6.0').
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/nextflow)

`--samples_tsv_file`
* A TSV (tab-separated values) file with the following headers:

| Header       | Description                   |
|--------------|-------------------------------|
| sample_id    | Sample ID.                    |
| fastq_file   | Full path to `fastq.gz` file. |

`--reference_genome_fasta_file`
* Reference genome FASTA (`.fasta` or `.fasta.gz`) file.

`--deepvariant_input_path`
* Parent path to input `fastq.gz` and `fasta` files.

`--deepvariant_output_path`
* Parent path to output directory.

`--deepvariant_containerization`
* Either `singularity` or `docker`.
