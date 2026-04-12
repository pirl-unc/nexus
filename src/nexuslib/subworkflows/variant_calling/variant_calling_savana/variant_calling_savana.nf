#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }                       from '../../../tools/samtools'
include { runSavanaRun }                           from '../../../tools/savana'
include { runSavanaClassify }                      from '../../../tools/savana'
include { decompressFile as decompressFasta }      from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir = ''
params.reference_genome_fasta_file      = ''
params.custom_params_file               = ''
params.contigs_txt_file                 = ''

// Optional arguments
params.params_savana_run                = '--length 30 --mapq 20 --min_support 3'
params.params_savana_classify           = ''

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_SAVANA {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)]
        reference_genome_fasta_file
        contigs_txt_file
        custom_params_file
        params_savana_run
        params_savana_classify
        output_dir

    main:
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidx(decompressFasta.out.f)
        fasta_file                  = runSamtoolsFaidx.out.fasta
        fasta_fai_file              = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file              = runSamtoolsFaidx.out.gzi_file

        runSavanaRun(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            contigs_txt_file,
            params_savana_run,
            output_dir
        )

        runSavanaClassify(
            runSavanaRun.out.f,
            custom_params_file,
            params_savana_classify,
            output_dir
        )

    emit:
        runSavanaClassify.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =======================================================================================
             Identify somatic structural variants in long-read DNA sequencing BAM files using Savana
             =======================================================================================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run Savana.

        usage: nexus run --nf-workflow variant_calling_savana.nf [required] [optional] [--help]

        required arguments:
            -c                                  :   Nextflow .config file.
            -w                                  :   Nextflow work directory path.
            --samples_tsv_file                  :   TSV file with the following columns:
                                                    'sample_id',
                                                    'tumor_bam_file',
                                                    'tumor_bam_bai_file',
                                                    'normal_bam_file',
                                                    'normal_bam_bai_file'
            --output_dir                        :   Directory to which output files will be copied.
            --reference_genome_fasta_file       :   Reference genome FASTA file.
            --custom_params_file                :   Savana classify --custom_params file.
            --contigs_txt_file                  :   TXT file with each contig name in a new line.

        optional arguments:
            --params_savana_run                 :   Savana run parameters (default: '"--length 30 --mapq 20 --min_support 3"').
                                                    Note that the parameters need to be wrapped in quotes.
            --params_savana_classify            :   Savana classify parameters (default: '""').
                                                    Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_savana_run      = (params.params_savana_run == true) ? '' : params.params_savana_run
    def params_savana_classify = (params.params_savana_classify == true) ? '' : params.params_savana_classify

    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        custom_params_file                  :   ${params.custom_params_file}
        contigs_txt_file                    :   ${params.contigs_txt_file}
        params_savana_run                   :   ${params_savana_run}
        params_savana_classify              :   ${params_savana_classify}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.tumor_bam_file}",
            "${row.tumor_bam_bai_file}",
            "${row.normal_bam_file}",
            "${row.normal_bam_bai_file}") }
        .set { input_bam_files_ch }

    VARIANT_CALLING_SAVANA(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.contigs_txt_file,
        params.custom_params_file,
        params_savana_run,
        params_savana_classify,
        params.output_dir
    )
}
