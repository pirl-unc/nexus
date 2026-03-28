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

if (params.params_vep == true) {
    params_vep = ''
} else {
    params_vep = params.params_vep
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
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
} else {
    log.info"""\
        samples_tsv_file        :   ${params.samples_tsv_file}
        vep_dir                 :   ${params.vep_dir}
        output_dir              :   ${params.output_dir}
        params_vep              :   ${params_vep}
    """.stripIndent()
}

// ------------------------------------------------------------
// Step 4. Set channels
// ------------------------------------------------------------
Channel
    .fromPath( params.samples_tsv_file )
    .splitCsv( header: true, sep: '\t' )
    .map { row -> tuple(
        "${row.sample_id}",
        "${row.vcf_file}") }
    .set { input_vcf_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
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
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_ANNOTATION_VEP(
        input_vcf_files_ch,
        params.vep_dir,
        params_vep,
        params.output_dir
    )
}
