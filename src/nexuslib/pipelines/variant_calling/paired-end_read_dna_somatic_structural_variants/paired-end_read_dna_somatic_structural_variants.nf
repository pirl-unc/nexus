#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runDelly2TumorNormal } from '../../modules/delly2'
include { runLumpyExpressTumorNormal } from '../../modules/lumpy'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.samples_tsv_file = ''
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.tools_list = 'delly2,lumpyexpress'
params.delly2 = 'delly'
params.delly2_params = '--exclude /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/delly2/human.hg38.excl.tsv --map-qual 20 '
params.bcftools = 'bcftools'
params.python2 = 'python2'
params.lumpyexpress = 'lumpyexpress'
params.lumpyexpress_config_file = 'lumpyexpress.config'
params.lumpy_extract_split_reads_script_file = '/datastore/lbcfs/collaborations/pirl/share/apps/lumpy/lumpy-sv-0.3.1/scripts/extractSplitReads_BwaMem'
params.samtools = 'samtools'
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         =========================================================================================
         Identify somatic structural variants using (Illumina) paired-end DNA sequencing bam files
         =========================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run delly2 (tumor and normal mode).
        2. Run lumpyexpress (tumor and normal mode).

    usage: nexus run --nf-workflow paired_end_dna_somatic_structural_variants.nf [required] [optional] [--help]

    required arguments:
        -c                                          :   Nextflow .config file.
        -w                                          :   Nextflow work directory path.
        --samples_tsv_file                          :   TSV file with the following columns:
                                                        'sample_id',
                                                        'tumor_bam_file',
                                                        'tumor_bam_bai_file',
                                                        'normal_bam_file',
                                                        'normal_bam_bai_file',
                                                        'tumor_sample_id',
                                                        'normal_sample_id'
        --output_dir                                :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file               :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --tools_list                                :   Tools to run (default: 'delly2,lumpyexpress').
        --delly2                                    :   delly2 path (default:
                                                        /datastore/lbcfs/collaborations/pirl/share/apps/delly2/delly_v1.1.8_linux_x86_64bit).
        --delly2_params                             :   delly2 'call' parameters (default:
                                                        '"--exclude /datastore/lbcfs/collaborations/pirl/seqdata/tool-resources/delly2/human.hg38.excl.tsv --map-qual 20 "').
                                                        Note that the parameters need to be wrapped in quotes
                                                        and a space at the end of the string is necessary.
        --bcftools                                  :   bcftools path (default: bcftools).
        --python2                                   :   python2 path (default: python2).
        --lumpyexpress                              :   lumpyexpress path (default:
                                                        /datastore/lbcfs/collaborations/pirl/share/apps/lumpy/lumpy-sv-0.3.1/scripts/lumpyexpress).
        --lumpyexpress_config_file                  :   lumpyexpress config file (default:
                                                        /datastore/lbcfs/collaborations/pirl/share/apps/lumpy/lumpy-sv-0.3.1/scripts/lumpyexpress.config).
        --lumpy_extract_split_reads_script_file     :   lumpy 'extractSplitReads_BwaMem' file (default:
                                                        /datastore/lbcfs/collaborations/pirl/share/apps/lumpy/lumpy-sv-0.3.1/scripts/extractSplitReads_BwaMem).
        --samtools                                  :   samtools path (default: samtools).
        --delete_work_dir                           :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                            :   ${params.samples_tsv_file}
        output_dir                                  :   ${params.output_dir}
        reference_genome_fasta_file                 :   ${params.reference_genome_fasta_file}
        tools_list                                  :   ${params.tools_list}
        delly2                                      :   ${params.delly2}
        delly2_params                               :   ${params.delly2_params}
        bcftools                                    :   ${params.bcftools}
        python2                                     :   ${params.python2}
        lumpyexpress                                :   ${params.lumpyexpress}
        lumpyexpress_config_file                    :   ${params.lumpyexpress_config_file}
        lumpy_extract_split_reads_script_file       :   ${params.lumpy_extract_split_reads_script_file}
        samtools                                    :   ${params.samtools}
        delete_work_dir                             :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Set channels
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.tumor_bam_file}",
        "${row.tumor_bam_bai_file}",
        "${row.normal_bam_file}",
        "${row.normal_bam_bai_file}",
        "${row.tumor_sample_id}",
        "${row.normal_sample_id}") }
    .set { input_bam_files_ch }

// Step 5. Workflow
workflow PAIRED_END_DNA_SOMATIC_STRUCTURAL_VARIANTS {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file), val(tumor_sample_id), val(normal_sample_id)]
        output_dir
        reference_genome_fasta_file
        tools_list
        delly2
        delly2_params
        bcftools
        python2
        lumpyexpress
        lumpyexpress_config_file
        lumpy_extract_split_reads_script_file
        samtools

    main:
        tools = tools_list?.split(',') as List

        input_bam_files_ch
            .map{ [it[0], it[1], it[2], it[3], it[4], it[5], it[6]] }
            .set{ run_delly2_input_ch }
        input_bam_files_ch
            .map{ [it[0], it[1], it[2], it[3], it[4], it[5], it[6]] }
            .set{ run_lumpy_input_ch }

        // delly2
        if (tools.contains('delly2')) {
            runDelly2TumorNormal(
                run_delly2_input_ch,
                delly2,
                delly2_params,
                bcftools,
                reference_genome_fasta_file,
                output_dir
            )
        }

        // lumpyexpress
        if (tools.contains('lumpyexpress')) {
            runLumpyExpressTumorNormal(
                run_lumpy_input_ch,
                python2,
                lumpyexpress,
                lumpyexpress_config_file,
                lumpy_extract_split_reads_script_file,
                samtools,
                output_dir
            )
        }
}

workflow {
    PAIRED_END_DNA_SOMATIC_STRUCTURAL_VARIANTS(
        input_bam_files_ch,
        params.output_dir,
        params.reference_genome_fasta_file,
        params.tools_list,
        params.delly2,
        params.delly2_params,
        params.bcftools,
        params.python2,
        params.lumpyexpress,
        params.lumpyexpress_config_file,
        params.lumpy_extract_split_reads_script_file,
        params.samtools
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
