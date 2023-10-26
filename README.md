# Available Pipelines

| Pipeline                            | Description                              |
|:------------------------------------|:-----------------------------------------|
| [Alignment](#alignment)             | Align sequencing reads |
| [Quantification](#quantification)   | Quantify analytes |
| [Variant Calling](#variant-calling) | Identify variants |

## Alignment
### [long_read_alignment_minimap2.nf](/pipelines/alignment/long_read_minimap2/)
Aligns long (Pacific Biosciences or Oxford Nanopore) DNA or RNA `FASTQ` files to a reference genome using `minimap2`.

### [long_read_rna_alignment_ultra.nf](/pipelines/alignment/long_read_rna_ultra/)
Aligns long (Pacific Biosciences or Oxford Nanopore) RNA `FASTQ` files  to a reference genome using `ultra`.

### [paired_end_read_dna_alignment_bwa-mem2.nf](/pipelines/alignment/paired_end_read_dna_bwa-mem2/)
Aligns paired-end (Illumina) DNA `FASTQ` files to a reference genome using `bwa-mem2`.

<hr/>

## Quantification
### [paired_end_read_rna_quantification_salmon_mapping_mode.nf](/pipelines/quantification/paired_end_read_rna/salmon_mapping_mode/)
Quantifies paired-end (Illumina) RNA `FASTQ` files using `salmon` (mapping mode).

<hr/>

## Variant Calling
### [pacbio_dna_small_variants.nf](/pipelines/variant_calling/long_read_dna_small_variants/pacbio/)
Identifies small DNA variants (SNVs and INDELs) using `deepvariant` in PacBio long-read DNA sequencing `BAM` files.

### [pacbio_dna_structural_variants.nf](/pipelines/variant_calling/long_read_dna_structural_variants/pacbio/)
Identifies structural DNA variants using `sniffles2`, `pbsv`, `svim`, and `cutesv` in PacBio long-read DNA sequencing `BAM` files.

### [pacbio_dna_somatic_structural_variants.nf](/pipelines/variant_calling/long_read_dna_structural_variants/pacbio/)
Identifies somatic structural DNA variants using `savana` in PacBio long-read DNA sequencing `BAM` files.

### [paired_end_human_dna_somatic_small_variants.nf](/pipelines/variant_calling/paired_end_dna_small_variants/somatic/human/)
Identifies small DNA somatic variants (SNVs and INDELs) using `gatk4-mutect2` and `strelka2` in paired-end Illumina human DNA sequencing `BAM` files.

### [paired_end_dna_somatic_structural_variants.nf](/pipelines/variant_calling/paired_end_dna_structural_variants/somatic/)
Identifies structural DNA somatic variants using `delly2` and `lumpyexpress` in paired-end Illumina DNA sequencing `BAM` files.

<hr/>
