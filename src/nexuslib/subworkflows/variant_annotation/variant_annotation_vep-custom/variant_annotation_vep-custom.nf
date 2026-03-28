#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }                from '../../../tools/samtools'
include { runVEPCustom }                    from '../../../tools/vep'
include { prepareGencodeGtfFileForVEP }     from '../../../tools/vep'

// ------------------------------------------------------------
// Step 2. Input parameters
// ------------------------------------------------------------
params.help                             = ''

// Required arguments
params.samples_tsv_file                 = ''
params.vep_dir                          = ''
params.output_dir                       = ''
params.reference_genome_fasta_file      = ''
params.reference_genes_gtf_file         = ''

// Optional arguments
params.reference_genes_gtf_source       = 'GENCODE'
params.params_vep                       = '--species homo_sapiens --offline --cache --assembly GRCh38'

if (params.params_vep == true) {
    params_vep = ''
} else {
    params_vep = params.params_vep
}

// ------------------------------------------------------------
// Step 3. Print inputs and help
// ------------------------------------------------------------
log.info """\
         ============================================================
         Annotate variants using VEP using custom FASTA and GTF files
         ============================================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run VEP.

    usage: nexus run --nf-workflow variant_annotation_vep-custom.nf [required] [optional] [--help]

    required arguments:
        -c                                  :   Nextflow .config file.
        -w                                  :   Nextflow work directory path.
        --samples_tsv_file                  :   TSV file with the following columns: 'sample_id', 'vcf_file'.
        --vep_dir                           :   VEP cache directory.
        --output_dir                        :   Directory to which output files will be copied.
        --reference_genome_fasta_file       :   Reference genome FASTA file.
        --reference_genes_gtf_file          :   Reference genes GTF file.

    optional arguments:
        --reference_genes_gtf_source        :   Reference genes GTF file source (default: GENCODE).
        --params_vep                        :   VEP parameters (default: '"--species homo_sapiens --database --offline --cache --assembly GRCh38"').
                                                Note that the parameters need to be wrapped in quotes.
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        samples_tsv_file                :   ${params.samples_tsv_file}
        vep_dir                         :   ${params.vep_dir}
        output_dir                      :   ${params.output_dir}
        reference_genome_fasta_file     :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file        :   ${params.reference_genes_gtf_file}
        reference_genes_gtf_source      :   ${params.reference_genes_gtf_source}
        params_vep                      :   ${params_vep}
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
workflow VARIANT_ANNOTATION_VEP_CUSTOM {
    take:
        input_vcf_files_ch             // channel: [val(sample_id), path(vcf_file)]
        vep_dir
        reference_genome_fasta_file
        reference_genes_gtf_file
        reference_genes_gtf_source
        params_vep
        output_dir

    main:
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file          = runSamtoolsFaidx.out.fasta
        fasta_fa_file       = runSamtoolsFaidx.out.fai_file

        prepareGencodeGtfFileForVEP(reference_genes_gtf_file)
        gtf_gz_file         = prepareGencodeGtfFileForVEP.out.gtf_file
        gtf_gz_tbi_file     = prepareGencodeGtfFileForVEP.out.tbi_file

        runVEPCustom(
            input_vcf_files_ch,
            vep_dir,
            fasta_file,
            fasta_fa_file,
            gtf_gz_file,
            gtf_gz_tbi_file,
            reference_genes_gtf_source,
            params_vep,
            output_dir
        )

    emit:
        runVEPCustom.out.f
}

// ------------------------------------------------------------
// Step 6. Entry workflow
// ------------------------------------------------------------
workflow {
    VARIANT_ANNOTATION_VEP_CUSTOM(
        input_vcf_files_ch,
        params.vep_dir,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params.reference_genes_gtf_source,
        params_vep,
        params.output_dir
    )
}
