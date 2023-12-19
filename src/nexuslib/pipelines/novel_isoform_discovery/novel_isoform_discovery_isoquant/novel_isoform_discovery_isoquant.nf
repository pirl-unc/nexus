//#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// Step 1. Import Nextflow modules
include { runIsoquant } from '../../modules/isoquant'

// Step 2. Input arguments
params.help = ''
// Required arguments
params.output_dir = ''
// Optional arguments
params.reference_genome_fasta_file = '/datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa'
params.isoquant = 'isoquant.py'
params.isoquant_params = '--data_type pacbio_ccs --sqanti_output --high_memory --complete_genedb --genedb /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf '
params.delete_work_dir = false

// Step 3. Print inputs and help
log.info """\
         ======================================
         Novel isoform discovery using Isoquant
         ======================================
         """.stripIndent()

if (params.help) {
    log.info"""\
    workflow:
        1. Run isoquant.

    usage: nexus run --nf-workflow novel_isoform_discovery_isoquant.nf [required] [optional] [--help]

    required arguments:
        -c                              :   Nextflow .config file.
        -w                              :   Nextflow work directory path.
        --output_dir                    :   Directory to which output files will be copied.

    optional arguments:
        --reference_genome_fasta_file   :   Reference genome FASTA file (default: /datastore/lbcfs/collaborations/pirl/seqdata/references/hg38.fa).
        --isoquant                      :   isoquant path (default: isoquant.py).
        --isoquant_params               :   isoquant parameters (default: '"--data_type pacbio_ccs --sqanti_output --high_memory --complete_genedb --genedb /datastore/lbcfs/collaborations/pirl/seqdata/references/gencode.v41.annotation.gtf "').
                                            Note that the parameters need to be wrapped in quotes
                                            and a space at the end of the string is necessary.
        --delete_work_dir               :   Delete work directory (default: false).
    """.stripIndent()
    exit 0
} else {
    log.info"""\
        output_dir                      :   ${params.output_dir}
        reference_genome_fasta_file     :   ${params.reference_genome_fasta_file}
        isoquant                        :   ${params.isoquant}
        isoquant_params                 :   ${params.isoquant_params}
        delete_work_dir                 :   ${params.delete_work_dir}
    """.stripIndent()
}

// Step 4. Workflow
workflow NOVEL_ISOFORM_DISCOVERY_ISOQUANT {
    take:
        reference_genome_fasta_file
        isoquant
        isoquant_params
        output_dir
    main:
        runIsoquant(
            reference_genome_fasta_file,
            isoquant,
            isoquant_params,
            output_dir
        )
    emit:
        runIsoquant.out.f
}

workflow {
    NOVEL_ISOFORM_DISCOVERY_ISOQUANT(
        params.reference_genome_fasta_file,
        params.isoquant,
        params.isoquant_params,
        params.output_dir
    )
}

workflow.onComplete {
    if ( params.delete_work_dir == true || params.delete_work_dir == 1 ) {
        workflow.workDir.deleteDir()
    }
}
