#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runVEP }  from '../../../tools/vep'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                 = ''

// Required arguments
params.samples_tsv_file     = ''
params.vep_dir              = ''
params.output_dir           = ''

// Optional arguments
params.params_vep           = '--species homo_sapiens --database --offline --cache --assembly GRCh38'

// ------------------------------------------------------------
// Step 3. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_ANNOTATION_VEP {
    take:
        input_vcf_files_ch             // channel: [val(sample_id), path(vcf_file)]
        vep_dir
        params_vep
        output_dir

    main:
        runVEP(
            input_vcf_files_ch,
            vep_dir,
            params_vep,
            output_dir
        )

    emit:
        runVEP.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===========================
             Annotate variants using VEP
             ===========================
             """.stripIndent()

    if (params.help) {
        log.info"""\
        workflow:
            1. Run VEP.

        usage: nexus run --nf-workflow variant_annotation_vep.nf [required] [optional] [--help]

        required arguments:
            -c                      :   Nextflow .config file.
            -w                      :   Nextflow work directory path.
            --samples_tsv_file      :   TSV file with the following columns: 'sample_id', 'vcf_file'.
            --vep_dir               :   VEP cache directory.
            --output_dir            :   Directory to which output files will be copied.

        optional arguments:
            --params_vep            :   VEP parameters (default: '"--species homo_sapiens --database --offline --cache --assembly GRCh38"').
                                        Note that the parameters need to be wrapped in quotes.
        """.stripIndent()
        exit 0
    }

    def params_vep = (params.params_vep == true) ? '' : params.params_vep

    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        vep_dir                 :   ${params.vep_dir}
        output_dir              :   ${params.output_dir}
        params_vep              :   ${params_vep}
    """.stripIndent()

    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.vcf_file}") }
        .set { input_vcf_files_ch }

    VARIANT_ANNOTATION_VEP(
        input_vcf_files_ch,
        params.vep_dir,
        params_vep,
        params.output_dir
    )
}
