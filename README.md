# Available Pipelines

| Pipeline                            | Description                              |
|:------------------------------------|:-----------------------------------------|
| [Alignment](#alignment)             | Align sequencing reads |
| [Quantification](#quantification)   | Quantify analytes |
| [Variant Calling](#variant-calling) | Identify variants |

## Alignment
### [long_read_alignment_minimap2.nf](/workflows/pipelines/alignment/long_read_alignment_minimap2/)
Aligns long DNA or RNA reads (Pacific Biosciences or Oxford Nanopore) to a reference genome using `minimap2`.

### [long_read_rna_alignment_ultra.nf](/workflows/pipelines/alignment/long_read_rna_alignment_ultra/)
Aligns long RNA reads (Pacific Biosciences or Oxford Nanopore) to a reference genome using `ultra`.

### [paired_end_read_dna_alignment_bwa-mem2.nf](/workflows/pipelines/alignment/paired_end_read_dna_alignment_bwa-mem2/)
Aligns paired-end DNA reads (Illumina) to a reference genome using `bwa-mem2`.

<hr/>

## Quantification
### [paired_end_read_rna_quantification_salmon_mapping_mode.nf](/workflows/pipelines/quantification/paired_end_read_rna/salmon_mapping_mode/)
Quantifies RNA using `salmon` (mapping mode).

<hr/>

## Variant Calling

<hr/>
