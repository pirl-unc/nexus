#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runSniffles2 } from '../../modules/sniffles2'
include { runPbsv } from '../../modules/pbsv'
include { runSvimAlignmentMode } from '../../modules/svim'
include { runCuteSV } from '../../modules/cutesv'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.tools_list = 'sniffles2,pbsv,svim,cutesv'
params.sniffles2 = 'sniffles'
params.sniffles2_params = '--minsupport 3 --minsvlen 30 --mapq 20 --output-rnames '
params.pbsv = 'pbsv'
params.pbsv_discover_params = '--ccs --min-gap-comp-id-perc 97 --min-mapq 20 '
params.pbsv_call_params = '--ccs --call-min-reads-per-strand-all-samples 0 --call-min-read-perc-one-sample 10 --call-min-reads-all-samples 3 --call-min-reads-one-sample 3 '
params.svim = 'svim'
params.svim_params = '--min_mapq 20 --min_sv_size 30 --insertion_sequences --read_names --zmws '
params.cutesv = 'cuteSV'
params.cutesv_params = '--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_support 3 --min_mapq 20 --min_size 30 --max_size -1 --report_readid --genotype '
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
        ==================================================================================
        Identify structural DNA variants using long-read whole-genome sequencing BAM files
        ==================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
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
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        output_dir                      :   ${params.output_dir}
        reference_genome_fasta_file     :   ${params.reference_genome_fasta_file}
        tools_list                      :   ${params.tools_list}
        sniffles2                       :   ${params.sniffles2}
        sniffles2_params                :   ${params.sniffles2_params}
        pbsv                            :   ${params.pbsv}
        pbsv_discover_params            :   ${params.pbsv_discover_params}
        pbsv_call_params                :   ${params.pbsv_call_params}
        svim                            :   ${params.svim}
        svim_params                     :   ${params.svim_params}
        cutesv                          :   ${params.cutesv}
        cutesv_params                   :   ${params.cutesv_params}
        delete_work_dir                 :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.bam_file}",
        "${row.bam_bai_file}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow LONG_READ_DNA_STRUCTURAL_VARIANTS {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(bam_file), path(bam_bai_file)]
        reference_genome_fasta_file
        tools_list
        sniffles2
        sniffles2_params
        pbsv
        pbsv_discover_params
        pbsv_call_params
        svim
        svim_params
        cutesv
        cutesv_params
        output_dir

    main:
        tools = tools_list?.split(',') as List

        run_sniffles2_input_ch = input_bam_files_ch
        run_pbsv_input_ch = input_bam_files_ch
        run_svim_input_ch = input_bam_files_ch
        run_cutesv_input_ch = input_bam_files_ch

        // Sniffles2
        if (tools.contains('sniffles2')) {
            runSniffles2(
                run_sniffles2_input_ch,
                reference_genome_fasta_file,
                sniffles2,
                sniffles2_params,
                output_dir
            )
            runSniffles2.out.f.set{ run_sniffles2_output_ch }
        }

        // pbsv
        if (tools.contains('pbsv')) {
            runPbsv(
                run_pbsv_input_ch,
                reference_genome_fasta_file,
                pbsv,
                pbsv_discover_params,
                pbsv_call_params,
                output_dir
            )
            runPbsv.out.f.set{ run_pbsv_output_ch }
        }

        // SVIM
        if (tools.contains('svim')) {
            runSvimAlignmentMode(
                run_svim_input_ch,
                reference_genome_fasta_file,
                svim,
                svim_params,
                output_dir
            )
            runSvimAlignmentMode.out.f.set{ run_svim_output_ch }
        }

        // cuteSV
        if (tools.contains('cutesv')) {
            runCuteSV(
                run_cutesv_input_ch,
                reference_genome_fasta_file,
                cutesv,
                cutesv_params,
                output_dir
            )
            runCuteSV.out.f.set{ run_cutesv_output_ch }
        }
}

workflow {
    LONG_READ_DNA_STRUCTURAL_VARIANTS(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.tools_list,
        params.sniffles2,
        params.sniffles2_params,
        params.pbsv,
        params.pbsv_discover_params,
        params.pbsv_call_params,
        params.svim,
        params.svim_params,
        params.cutesv,
        params.cutesv_params,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
