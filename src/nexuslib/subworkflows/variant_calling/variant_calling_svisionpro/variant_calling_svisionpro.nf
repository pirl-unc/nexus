#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }        from '../../../tools/samtools'
include { runSVisionPro }           from '../../../tools/svisionpro'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help = ''

// Required arguments
params.samples_tsv_file                 = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.svisionpro_model_file            = ''

// Optional arguments
params.params_svisionpro                = '--detect_mode somatic --preset hifi --min_supp 3 --min_mapq 20 --min_sv_size 30 --max_sv_size 1000000 --device cpu --img_size 256'
params.params_svisionpro_extract        = '--extract somatic --min_supp 3'

if (params.params_svisionpro == true) {
    params_svisionpro = ''
} else {
    params_svisionpro = params.params_svisionpro
}

if (params.params_svisionpro_extract == true) {
    params_svisionpro_extract = ''
} else {
    params_svisionpro_extract = params.params_svisionpro_extract
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ===========================================================================================
         Identify somatic structural variants in long-read DNA sequencing BAM files using SVisionPro
         ===========================================================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run SVision-pro.

    usage: nexus run --nf-workflow variant_calling_svisionpro.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'tumor_bam_file', 'tumor_bam_bai_file', 'normal_bam_file', 'normal_bam_bai_file'.
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.
        --svisionpro_model_file             :   SVision-pro model file.

    optional arguments:
        --params_svisionpro                 :   SVision-pro parameters (default: '"--detect_mode somatic --preset hifi --min_supp 3 --min_mapq 20 --min_sv_size 30 --max_sv_size 1000000 --device cpu --img_size 256"').
                                                Note that the parameters need to be wrapped in quotes.
        --params_svisionpro_extract         :   SVision-pro extract_op.py parameters (default: '"--extract somatic --min_supp 3"').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                    :   ${params.samples_tsv_file}
        output_dir                          :   ${params.output_dir}
        reference_genome_fasta_file         :   ${params.reference_genome_fasta_file}
        svisionpro_model_file               :   ${params.svisionpro_model_file}
        params_svisionpro                   :   ${params_svisionpro}
        params_svisionpro_extract           :   ${params_svisionpro_extract}
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
        "${row.tumor_bam_file}",
        "${row.tumor_bam_bai_file}",
        "${row.normal_bam_file}",
        "${row.normal_bam_bai_file}") }
    .set { input_bam_files_ch }

// ------------------------------------------------------------
// Step 5. Sub-workflows
// ------------------------------------------------------------
workflow VARIANT_CALLING_SVISIONPRO {
    take:
        input_bam_files_ch             // channel: [val(sample_id), path(tumor_bam_file), path(tumor_bam_bai_file), path(normal_bam_file), path(normal_bam_bai_file)]
        reference_genome_fasta_file
        svisionpro_model_file
        params_svisionpro
        params_svisionpro_extract
        output_dir

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file          = runSamtoolsFaidx.out.fasta
        fasta_fai_file      = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file      = runSamtoolsFaidx.out.gzi_file

        runSVisionPro(
            input_bam_files_ch,
            fasta_file,
            fasta_fai_file,
            svisionpro_model_file,
            params_svisionpro,
            params_svisionpro_extract,
            output_dir
        )

    emit:
        runSVisionPro.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_CALLING_SVISIONPRO(
        input_bam_files_ch,
        params.reference_genome_fasta_file,
        params.svisionpro_model_file,
        params_svisionpro,
        params_svisionpro_extract,
        params.output_dir
    )
}
