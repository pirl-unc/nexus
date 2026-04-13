#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }              from '../../../tools/samtools'
include { runSamtoolsFaidxFasta }         from '../../../tools/samtools'
include { prepareGencodeGtfFileForVEP }   from '../../../tools/vep'
include { runVEPCustom }                  from '../../../tools/vep'
include { runFilterVEP }                  from '../../../tools/vep'
include { runMoPepGenParseVEP }           from '../../../tools/mopepgen'
include { runMoPepGenParseREDItools2 }    from '../../../tools/mopepgen'
include { runMoPepGenParseArriba }        from '../../../tools/mopepgen'
include { runMoPepGenParseRMATS }         from '../../../tools/mopepgen'
include { runMoPepGenParseCircExplorer2 } from '../../../tools/mopepgen'
include { runMoPepGenCallVariant }        from '../../../tools/mopepgen'
include { decompressFile as decompressFasta }      from '../../../tools/utils'
include { decompressFile as decompressGtf }        from '../../../tools/utils'
include { decompressFile as decompressProteome }   from '../../../tools/utils'

// ------------------------------------------------------------
// Step 2. Helper functions
// ------------------------------------------------------------
def isValidFile(filePath) {
    if (filePath == null || filePath == '' || filePath == 'NA' || filePath == 'na') return false
    return file(filePath).exists()
}

// ------------------------------------------------------------
// Step 3. Workflow
// ------------------------------------------------------------
workflow PEPTIDE_PREDICTION_MOPEPGEN {
    take:
        vep_input_ch                // channel: [val(sample_id), path(vcf_file)]
        caller_metadata_ch          // channel: [val(sample_id), val(vcf_basename), val(caller)]
        reditools2_ch               // channel: [val(sample_id), path(tsv_file)]
        arriba_ch                   // channel: [val(sample_id), path(tsv_file)]
        rmats_ch                    // channel: [val(sample_id), path(rmats_dir)]
        circexplorer2_ch            // channel: [val(sample_id), path(bed_file)]
        reference_genome_fasta_file
        reference_genes_gtf_file
        reference_proteome_fasta_file
        output_dir
        cfg_vep
        cfg_mopepgen

    main:
        // ---- Prepare reference files ----
        // Uncompressed FASTA (for moPepGen processes)
        decompressFasta(reference_genome_fasta_file)
        uncompressed_fasta = decompressFasta.out.f

        // Uncompressed GTF (for moPepGen processes)
        decompressGtf(reference_genes_gtf_file)
        uncompressed_gtf = decompressGtf.out.f

        // Uncompressed proteome FASTA (for moPepGen callVariant)
        decompressProteome(reference_proteome_fasta_file)
        uncompressed_proteome = decompressProteome.out.f

        // FASTA + index (for VEP)
        runSamtoolsFaidxFasta(reference_genome_fasta_file)
        genome_fasta     = runSamtoolsFaidxFasta.out.fasta
        genome_fasta_fai = runSamtoolsFaidxFasta.out.fai_file

        // Bgzipped + indexed GTF (for VEP)
        prepareGencodeGtfFileForVEP(reference_genes_gtf_file)
        vep_gtf_gz  = prepareGencodeGtfFileForVEP.out.gtf_file
        vep_gtf_tbi = prepareGencodeGtfFileForVEP.out.tbi_file

        // ---- VEP annotation (all VCFs flow through one invocation) ----
        // Uses real sample_id so publishDir outputs to output_dir/sample_id/
        runVEPCustom(
            vep_input_ch,
            cfg_vep.dir,
            genome_fasta,
            genome_fasta_fai,
            vep_gtf_gz,
            vep_gtf_tbi,
            cfg_vep.reference_source,
            cfg_vep.extra_args ?: '',
            output_dir
        )

        // ---- Filter VEP results ----
        runFilterVEP(
            runVEPCustom.out.f,
            cfg_vep.filter_extra_args ?: '',
            output_dir
        )

        // ---- Rejoin caller metadata with filtered VEP results ----
        // filterVEP output: tuple(sample_id, filtered_tsv)
        // filtered_tsv basename follows pattern: {original_vcf_basename}_vep-custom_filtered
        // We extract the original vcf_basename to join with caller_metadata_ch
        filtered_with_basename_ch = runFilterVEP.out.f.map { sid, filtered_tsv ->
            // Strip _vep-custom_filtered suffix to recover original VCF basename
            def basename = filtered_tsv.baseName
                .replaceAll('_vep-custom_filtered$', '')
                .replaceAll('_vep-custom$', '')
            tuple(sid, basename, filtered_tsv)
        }

        // Join on (sample_id, vcf_basename) to recover caller name
        // caller_metadata_ch: tuple(sample_id, vcf_basename, caller)
        // filtered_with_basename_ch: tuple(sample_id, vcf_basename, filtered_tsv)
        parsevep_joined_ch = filtered_with_basename_ch
            .combine(caller_metadata_ch, by: [0, 1])
            // Result: tuple(sample_id, vcf_basename, filtered_tsv, caller)

        parsevep_input_ch = parsevep_joined_ch.map { sid, basename, tsv, caller ->
            tuple(sid, tsv)
        }
        parsevep_source_ch = parsevep_joined_ch.map { sid, basename, tsv, caller ->
            caller
        }

        runMoPepGenParseVEP(
            parsevep_input_ch,
            uncompressed_fasta,
            uncompressed_gtf,
            cfg_vep.reference_source,
            parsevep_source_ch,
            cfg_mopepgen.parsevep_extra_args ?: '',
            output_dir
        )

        // ---- Direct parsers ----
        runMoPepGenParseREDItools2(
            reditools2_ch,
            uncompressed_gtf,
            cfg_vep.reference_source,
            cfg_mopepgen.parsereditools_extra_args ?: '',
            output_dir
        )

        runMoPepGenParseArriba(
            arriba_ch,
            uncompressed_fasta,
            uncompressed_gtf,
            cfg_vep.reference_source,
            cfg_mopepgen.parsearriba_extra_args ?: '',
            output_dir
        )

        runMoPepGenParseRMATS(
            rmats_ch,
            uncompressed_fasta,
            uncompressed_gtf,
            cfg_vep.reference_source,
            cfg_mopepgen.parsermats_extra_args ?: '',
            output_dir
        )

        runMoPepGenParseCircExplorer2(
            circexplorer2_ch,
            uncompressed_gtf,
            cfg_vep.reference_source,
            cfg_mopepgen.parsecircexplorer2_extra_args ?: '',
            output_dir
        )

        // ---- Collect all GVF results and group by sample_id ----
        all_gvf_ch = runMoPepGenParseVEP.out.f
            .mix(runMoPepGenParseREDItools2.out.f)
            .mix(runMoPepGenParseArriba.out.f)
            .mix(runMoPepGenParseRMATS.out.f)
            .mix(runMoPepGenParseCircExplorer2.out.f)
            .groupTuple()   // tuple(sample_id, [gvf1, gvf2, ...])

        // ---- Call variants ----
        runMoPepGenCallVariant(
            all_gvf_ch,
            uncompressed_fasta,
            uncompressed_gtf,
            uncompressed_proteome,
            cfg_mopepgen.callvariant_extra_args ?: '',
            output_dir
        )

    emit:
        runMoPepGenCallVariant.out.f
}

// ------------------------------------------------------------
// Step 4. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             ===============================================
             Predict mutant peptide sequences using moPepGen
             ===============================================
             """.stripIndent()

    if (params.help) {
        log.info """\
        usage: nexus run --nf-workflow peptide_prediction_mopepgen-gencode.nf -params-file params.yaml [--help]

        All parameters are supplied via a params.yaml file. See params.yaml for
        full documentation and defaults.
        """.stripIndent()
        exit 0
    }

    if (!params.samples_tsv_file)              error "ERROR: samples_tsv_file is required."
    if (!params.output_dir)                    error "ERROR: output_dir is required."
    if (!params.reference_genome_fasta_file)   error "ERROR: reference_genome_fasta_file is required."
    if (!params.reference_genes_gtf_file)      error "ERROR: reference_genes_gtf_file is required."
    if (!params.reference_proteome_fasta_file) error "ERROR: reference_proteome_fasta_file is required."
    if (!params.vep.dir)                       error "ERROR: vep.dir is required."

    log.info """\
        samples_tsv_file             :   ${params.samples_tsv_file}
        reference_genome_fasta_file  :   ${params.reference_genome_fasta_file}
        reference_genes_gtf_file     :   ${params.reference_genes_gtf_file}
        reference_proteome_fasta_file:   ${params.reference_proteome_fasta_file}
        output_dir                   :   ${params.output_dir}
        vep.dir                      :   ${params.vep.dir}
        """.stripIndent()

    // Auto-discover VCF columns (any column ending in _vcf_file).
    // Emit tuple(sample_id, caller_name, vcf_file) for each valid VCF.
    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .flatMap { row ->
            def entries = []
            row.each { col, val ->
                if (col.endsWith('_vcf_file') && isValidFile(val)) {
                    def caller = col.replace('_vcf_file', '')
                    entries.add(tuple(row.sample_id, caller, val))
                }
            }
            return entries
        }
        .set { vep_input_with_caller_ch }

    // Split into: VEP input (sample_id, vcf_file) and caller metadata (sample_id, caller, vcf_basename)
    def vep_input_ch = vep_input_with_caller_ch.map { sid, caller, vcf -> tuple(sid, vcf) }
    // Track caller name keyed by (sample_id, vcf_basename) so we can rejoin after VEP/filterVEP
    def caller_metadata_ch = vep_input_with_caller_ch.map { sid, caller, vcf ->
        tuple(sid, file(vcf).baseName, caller)
    }

    // Parse direct-parser columns from the same TSV
    Channel
        .fromPath( params.samples_tsv_file )
        .splitCsv( header: true, sep: '\t' )
        .set { tsv_rows_ch }

    // Direct parser input channels (only emitted when file is valid)
    def reditools2_ch    = tsv_rows_ch.filter { isValidFile(it.reditools2_tsv_file) }
                                      .map { row -> tuple(row.sample_id, row.reditools2_tsv_file) }
    def arriba_ch        = tsv_rows_ch.filter { isValidFile(it.arriba_tsv_file) }
                                      .map { row -> tuple(row.sample_id, row.arriba_tsv_file) }
    def rmats_ch         = tsv_rows_ch.filter { isValidFile(it.rmats_output_dir) }
                                      .map { row -> tuple(row.sample_id, row.rmats_output_dir) }
    def circexplorer2_ch = tsv_rows_ch.filter { isValidFile(it.circexplorer2_bed_file) }
                                      .map { row -> tuple(row.sample_id, row.circexplorer2_bed_file) }

    PEPTIDE_PREDICTION_MOPEPGEN(
        vep_input_ch,
        caller_metadata_ch,
        reditools2_ch,
        arriba_ch,
        rmats_ch,
        circexplorer2_ch,
        params.reference_genome_fasta_file,
        params.reference_genes_gtf_file,
        params.reference_proteome_fasta_file,
        params.output_dir,
        params.vep,
        params.mopepgen
    )
}
