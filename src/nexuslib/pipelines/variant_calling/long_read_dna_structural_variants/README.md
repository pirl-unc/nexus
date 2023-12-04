## long_read_dna_structural_variants.nf

Identifies structural DNA variants with long-read DNA sequencing 
(BAM) files using 
[Sniffles2](https://github.com/fritzsedlazeck/Sniffles), 
[Pbsv](https://github.com/PacificBiosciences/pbsv), 
[SVIM](https://github.com/eldariont/svim), and 
[cuteSV](https://github.com/tjiangHIT/cuteSV).

### Inputs / Outputs

| I/O    | Description                  |
|:-------|:-----------------------------|
| Input  | `bam` file for each sample.  | 
| Output | `vcf` files for each sample. |

### Dependencies

* `sniffles2`
* `pbsv`
* `svim`
* `cutesv`

### Usage

```
workflow:
    1. Run Sniffles2.
    2. Run PBSV.
    3. Run SVIM.
    4. Run cuteSV.

usage: nexus run --nf-workflow long_read_dna_structural_variants.nf [required] [optional] [--help]

required arguments:
    -c                              :   Nextflow .config file.
    -w                              :   Nextflow work directory path.
    --samples_tsv_file              :   TSV file with the following columns:
                                        'sample_id', 'bam_file', 'bam_bai_file'.
    --output_dir                    :   Directory to which output files will be copied.

optional arguments:
    --reference_genome_fasta_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
    --tools_list                    :   Tools to run (default: 'sniffles2,pbsv,svim,cutesv').
    --sniffles2                     :   sniffles2 path (default: sniffles).
    --sniffles2_params              :   sniffles2 parameters (default: '"--minsupport 3 --minsvlen 30 --mapq 20 --output-rnames "').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
    --pbsv                          :   pbsv path (default: pbsv).
    --pbsv_discover_params          :   pbsv 'discover' parameters (default: '"--ccs --min-gap-comp-id-perc 97 --min-mapq 20 "').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
    --pbsv_call_params              :   pbsv 'call' parameters (default:
                                        '"--ccs
                                          --call-min-reads-per-strand-all-samples 0
                                          --call-min-read-perc-one-sample 10
                                          --call-min-reads-all-samples 3
                                          --call-min-reads-one-sample 3 "').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
    --svim                          :   svim path (default: svim).
    --svim_params                   :   svim parameters (default: '"--min_mapq 20 --min_sv_size 30 --insertion_sequences --read_names --zmws "').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
    --cutesv                        :   cutesv path (default: cuteSV).
    --cutesv_params                 :   cutesv parameters (default:
                                        '"--max_cluster_bias_INS 1000
                                          --diff_ratio_merging_INS 0.9
                                          --max_cluster_bias_DEL 1000
                                          --diff_ratio_merging_DEL 0.5
                                          --min_support 3
                                          --min_mapq 20
                                          --min_size 30
                                          --max_size -1
                                          --report_readid
                                          --genotype "').
                                        Note that the parameters need to be wrapped in quotes
                                        and a space at the end of the string is necessary.
    --delete_work_dir               :   Delete work directory (default: false).
```

### Parameters

`-c`
* Nextflow config file can be downloaded [here](https://github.com/pirl-unc/nexus/tree/main/configs)

`--sample_tsv_file`

| Header       | Description                 |
|--------------|-----------------------------|
| sample_id    | Sample ID                   |
| bam_file     | Full path to `bam` file     |
| bam_bai_file | Full path to `bam.bai` file |
