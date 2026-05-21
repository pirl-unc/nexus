#!/usr/bin/env nextflow

/*
 * Author: Jin Seok (Andy) Lee
 *
 * Predict mutant proteoforms from raw long-read FASTQ.GZ files using exacto
 * v0.4.6a1. Mirrors the reference pipeline:
 *   https://pirl-unc.github.io/exacto/pipelines/mutant-proteoform-prediction.html
 *
 *   Per sample:
 *     1. minimap2 (map-hifi)            tumor DNA fastq.gz  → tumor DNA BAM
 *     2. minimap2 (map-hifi)            normal DNA fastq.gz → normal DNA BAM
 *     3. exacto call-somatic-dna-vars   tumor + normal DNA BAMs → somatic TSV
 *     4. exacto annotate-vars           somatic TSV + GTF → annotated DNA TSV
 *     5. RNA-Bloom2                     tumor RNA fastq.gz → assembled transcripts
 *     6. nexus_filter_rnabloom2_transcripts → filtered transcripts FASTA
 *     7. minimap2 (splice:hq)           filtered transcripts FASTA → transcriptome BAM
 *     8. exacto remove-unspliced-rnas   transcriptome BAM + FASTA + GTF → spliced BAM + FASTA
 *     9. exacto call-rna-vars           spliced BAM → RNA variants + transcript structures
 *    10. exacto integrate-vars          annotated DNA TSV + RNA variants TSV → integrated TSV
 *    11. exacto translate-structs       transcript structures + RNA vars + integrated vars → primary structures
 *    12. exacto call-peptide-vars       primary structures + reference proteome → mutant peptide variants
 */

nextflow.enable.dsl=2

// ------------------------------------------------------------
// Step 1. Import Nextflow modules
// ------------------------------------------------------------
include { runSamtoolsFaidx }                                  from '../../../tools/samtools'
include { runSamtoolsFaidxFasta }                             from '../../../tools/samtools'
include { decompressFile as decompressFasta }                 from '../../../tools/utils'
include { runMinimap2SortedBam as runMinimap2TumorDNA }       from '../../../tools/minimap2'
include { runMinimap2SortedBam as runMinimap2NormalDNA }      from '../../../tools/minimap2'
include { runMinimap2SortedBam as runMinimap2TumorRNA }       from '../../../tools/minimap2'
include { runRnaBloom2LongRead }                              from '../../../tools/rnabloom2'
include { runNexusFilterRNABloom2Transcripts }                from '../../../tools/nexus'
include { runExactoRemoveUnsplicedRNAs }                      from '../../../tools/exacto'
include { runExactoCallRNAVars }                              from '../../../tools/exacto'
include { runExactoCallSomaticDNAVars }                       from '../../../tools/exacto'
include { runExactoAnnotateVars }                             from '../../../tools/exacto'
include { runExactoIntegrateVars }                            from '../../../tools/exacto'
include { runExactoTranslateStructs }                         from '../../../tools/exacto'
include { runExactoCallPeptideVars }                          from '../../../tools/exacto'

// ------------------------------------------------------------
// Step 2. Workflow
// ------------------------------------------------------------
workflow PEPTIDE_PREDICTION_EXACTO {
    take:
        // samples_ch: [val(sample_id), val(tumor_dna_fastq_file),
        //              val(normal_dna_fastq_file), val(tumor_rna_fastq_file)]
        samples_ch
        reference_genome_fasta_file
        reference_gene_annotation_file
        reference_gene_annotation_source
        reference_gene_annotation_assembly
        reference_gene_annotation_version
        reference_proteome_fasta_file
        strategy
        platform_tag
        platform_unit_tag
        library_tag
        output_dir
        cfg

    main:
        // ---- Index reference genome FASTA (bgzipped) for minimap2 ----
        runSamtoolsFaidx(reference_genome_fasta_file)
        fasta_file     = runSamtoolsFaidx.out.fasta
        fasta_fai_file = runSamtoolsFaidx.out.fai_file
        fasta_gzi_file = runSamtoolsFaidx.out.gzi_file

        // ---- Decompressed reference + index (required by samtools calmd) ----
        decompressFasta(reference_genome_fasta_file)
        runSamtoolsFaidxFasta(decompressFasta.out.f)
        uncompressed_fasta     = runSamtoolsFaidxFasta.out.fasta
        uncompressed_fasta_fai = runSamtoolsFaidxFasta.out.fai_file

        // ---- Split samples row into three per-fastq sub-channels.
        //      Each sub-channel uses a SUFFIXED sample_id (e.g.
        //      "<sid>_tumor_dna") so minimap2 / publishDir outputs from the
        //      three alignment jobs do not collide. Suffixes are stripped
        //      again before downstream exacto steps that need to be joined
        //      back to a single base sample_id. ----
        tumor_dna_ch = samples_ch.map { sid, tdna_fq, ndna_fq, trna_fq ->
            tuple("${sid}_tumor_dna", [file(tdna_fq)])
        }
        normal_dna_ch = samples_ch.map { sid, tdna_fq, ndna_fq, trna_fq ->
            tuple("${sid}_normal_dna", [file(ndna_fq)])
        }
        tumor_rna_ch = samples_ch.map { sid, tdna_fq, ndna_fq, trna_fq ->
            tuple("${sid}_tumor_rna", file(trna_fq))
        }

        // ---- Step 1. Align tumor DNA (map-hifi) ----
        runMinimap2TumorDNA(
            tumor_dna_ch,
            fasta_file,
            fasta_fai_file,
            uncompressed_fasta,
            uncompressed_fasta_fai,
            cfg.minimap2_dna_args,
            platform_tag,
            platform_unit_tag,
            library_tag
        )

        // ---- Step 2. Align normal DNA (map-hifi) ----
        runMinimap2NormalDNA(
            normal_dna_ch,
            fasta_file,
            fasta_fai_file,
            uncompressed_fasta,
            uncompressed_fasta_fai,
            cfg.minimap2_dna_args,
            platform_tag,
            platform_unit_tag,
            library_tag
        )

        // ---- Step 3. exacto call-somatic-dna-vars (tumor vs matched normal) ----
        tumor_dna_bam_ch = runMinimap2TumorDNA.out.f.map { suffixed_sid, bam, bai ->
            def base = suffixed_sid.replaceAll(/_tumor_dna$/, '')
            tuple(base, bam, bai)
        }
        normal_dna_bam_ch = runMinimap2NormalDNA.out.f.map { suffixed_sid, bam, bai ->
            def base = suffixed_sid.replaceAll(/_normal_dna$/, '')
            tuple(base, bam, bai)
        }
        somatic_input_ch = tumor_dna_bam_ch.join(normal_dna_bam_ch)
            .map { sid, tbam, tbai, nbam, nbai ->
                tuple(sid, tbam, tbai, nbam, nbai)
            }
        runExactoCallSomaticDNAVars(
            somatic_input_ch,
            fasta_file,
            fasta_fai_file,
            fasta_gzi_file,
            cfg.call_somatic_dna_vars_extra_args,
            output_dir
        )

        // ---- Step 4. exacto annotate-vars (somatic DNA variants) ----
        runExactoAnnotateVars(
            runExactoCallSomaticDNAVars.out.f,
            reference_gene_annotation_file,
            reference_gene_annotation_source,
            reference_gene_annotation_assembly,
            reference_gene_annotation_version,
            cfg.annotate_vars_extra_args,
            output_dir
        )

        // ---- Step 5. RNA-Bloom2 assembly of tumor RNA reads ----
        runRnaBloom2LongRead(
            tumor_rna_ch,
            cfg.rnabloom2_extra_args,
            output_dir
        )

        // ---- Step 6. nexus_filter_rnabloom2_transcripts ----
        runNexusFilterRNABloom2Transcripts(
            runRnaBloom2LongRead.out.f,
            cfg.filter_rnabloom2_extra_args,
            output_dir
        )

        // ---- Step 7. Align filtered transcripts FASTQ.GZ to genome (splice:hq) ----
        // Use the FASTQ.gz (with per-base quality strings) rather than the FASTA
        // so the downstream BAM consumed by call-rna-vars carries QUAL fields.
        rnabloom2_fastq_ch = runNexusFilterRNABloom2Transcripts.out.f.map {
            suffixed_sid, reads_tsv, transcripts_tsv, transcripts_fa, transcripts_fq ->
                tuple(suffixed_sid, [transcripts_fq])
        }
        runMinimap2TumorRNA(
            rnabloom2_fastq_ch,
            fasta_file,
            fasta_fai_file,
            uncompressed_fasta,
            uncompressed_fasta_fai,
            cfg.minimap2_rna_args,
            platform_tag,
            platform_unit_tag,
            library_tag
        )

        // ---- Step 8. exacto remove-unspliced-rnas ----
        // Re-pair the transcriptome BAM with the transcripts FASTA, then
        // strip the _tumor_rna suffix so downstream sample_ids match the
        // base sample_id used by DNA somatic / annotate steps.
        rna_bam_ch = runMinimap2TumorRNA.out.f      // [suffixed_sid, bam, bai]
        rna_fasta_only_ch = runNexusFilterRNABloom2Transcripts.out.f.map {
            suffixed_sid, reads_tsv, transcripts_tsv, transcripts_fa, transcripts_fq ->
                tuple(suffixed_sid, transcripts_fa)
        }
        rna_assembly_ch = rna_bam_ch
            .join(rna_fasta_only_ch)
            .map { suffixed_sid, bam, bai, transcripts_fa ->
                def base = suffixed_sid.replaceAll(/_tumor_rna$/, '')
                tuple(base, bam, bai, transcripts_fa)
            }
        runExactoRemoveUnsplicedRNAs(
            rna_assembly_ch,
            reference_gene_annotation_file,
            reference_gene_annotation_source,
            reference_gene_annotation_assembly,
            reference_gene_annotation_version,
            cfg.remove_unspliced_rnas_extra_args ?: '',
            output_dir
        )

        // ---- Step 9. exacto call-rna-vars (on spliced BAM) ----
        rna_var_input_ch = runExactoRemoveUnsplicedRNAs.out.f.map { sid, bam, bai, _fa ->
            tuple(sid, bam, bai)
        }
        runExactoCallRNAVars(
            rna_var_input_ch,
            fasta_file,
            fasta_fai_file,
            fasta_gzi_file,
            reference_gene_annotation_file,
            reference_gene_annotation_source,
            reference_gene_annotation_assembly,
            reference_gene_annotation_version,
            cfg.call_rna_vars_extra_args,
            output_dir
        )

        // ---- Step 10. exacto integrate-vars (DNA + RNA) ----
        rna_calls_only_ch = runExactoCallRNAVars.out.f.map {
            sid, exons, rfs, rtm, introns, trsr, transcripts, transcript_structs, rna_calls ->
                tuple(sid, rna_calls)
        }
        integrate_input_ch = runExactoAnnotateVars.out.f
            .join(rna_calls_only_ch)
            .map { sid, ann_dna_tsv, rna_tsv ->
                // Use sample_id for both dna_id / rna_id by default
                tuple(sid, sid, sid, ann_dna_tsv, rna_tsv)
            }
        runExactoIntegrateVars(
            integrate_input_ch,
            reference_gene_annotation_file,
            reference_gene_annotation_source,
            reference_gene_annotation_assembly,
            reference_gene_annotation_version,
            cfg.integrate_vars_extra_args,
            output_dir
        )

        // ---- Step 11. exacto translate-structs → primary structures ----
        // translate-structs consumes the *transcript_structures* TSV (per-row
        // exon/junction structure) -- NOT the transcripts.tsv (per-transcript
        // summary) which lives at index `transcripts` in the call-rna-vars output.
        transcripts_and_rna_calls_ch = runExactoCallRNAVars.out.f.map {
            sid, exons, rfs, rtm, introns, trsr, transcripts, transcript_structs, rna_calls ->
                tuple(sid, transcript_structs, rna_calls)
        }
        integrated_only_ch = runExactoIntegrateVars.out.f.map {
            sid, dna_id, rna_id, integ_tsv -> tuple(sid, integ_tsv)
        }
        translate_input_ch = transcripts_and_rna_calls_ch
            .join(integrated_only_ch)
            .map { sid, transcripts, rna_calls, integ_tsv ->
                tuple(sid, transcripts, rna_calls, integ_tsv)
            }
        runExactoTranslateStructs(
            translate_input_ch,
            strategy,
            cfg.translate_structs_extra_args,
            output_dir
        )

        // ---- Step 12. exacto call-peptide-vars (vs reference proteome) ----
        peptide_input_ch = runExactoTranslateStructs.out.f.map {
            sid, primary_tsv, _primary_fa -> tuple(sid, primary_tsv)
        }
        runExactoCallPeptideVars(
            peptide_input_ch,
            reference_proteome_fasta_file,
            cfg.call_peptide_vars_extra_args,
            output_dir
        )

    emit:
        runExactoCallPeptideVars.out.f
}

// ------------------------------------------------------------
// Step 3. Entry workflow (runs only when this file is the main script)
// ------------------------------------------------------------
workflow {
    log.info """\
             =======================================
             Predict mutant proteoforms using Exacto
             =======================================
             """.stripIndent()

    if (params.help) {
        log.info """\
        usage: nexus run --nf-workflow peptide_prediction_exacto.nf -params-file params.yaml [--help]

        All parameters are supplied via a params.yaml file. See params.yaml for
        full documentation and defaults.
        """.stripIndent()
        exit 0
    }

    if (!params.samples_tsv_file)               error "ERROR: samples_tsv_file is required."
    if (!params.output_dir)                     error "ERROR: output_dir is required."
    if (!params.reference_genome_fasta_file)    error "ERROR: reference_genome_fasta_file is required."
    if (!params.reference_gene_annotation_file) error "ERROR: reference_gene_annotation_file is required."
    if (!params.reference_proteome_fasta_file)  error "ERROR: reference_proteome_fasta_file is required."

    log.info """\
        samples_tsv_file                     :   ${params.samples_tsv_file}
        output_dir                           :   ${params.output_dir}
        reference_genome_fasta_file          :   ${params.reference_genome_fasta_file}
        reference_gene_annotation_file       :   ${params.reference_gene_annotation_file}
        reference_gene_annotation_source     :   ${params.reference_gene_annotation_source}
        reference_gene_annotation_assembly   :   ${params.reference_gene_annotation_assembly}
        reference_gene_annotation_version    :   ${params.reference_gene_annotation_version}
        reference_proteome_fasta_file        :   ${params.reference_proteome_fasta_file}
        strategy                             :   ${params.strategy}
        platform_tag                         :   ${params.platform_tag}
        platform_unit_tag                    :   ${params.platform_unit_tag}
        library_tag                          :   ${params.library_tag}
        """.stripIndent()

    // ---- Samples TSV columns ----
    //   sample_id
    //   tumor_dna_fastq_file       long-read tumor DNA fastq.gz
    //   normal_dna_fastq_file      long-read matched-normal DNA fastq.gz
    //   tumor_rna_fastq_file       long-read tumor RNA fastq.gz
    Channel.fromPath(params.samples_tsv_file)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(
            "${row.sample_id}",
            "${row.tumor_dna_fastq_file}",
            "${row.normal_dna_fastq_file}",
            "${row.tumor_rna_fastq_file}") }
        .set { samples_ch }

    PEPTIDE_PREDICTION_EXACTO(
        samples_ch,
        params.reference_genome_fasta_file,
        params.reference_gene_annotation_file,
        params.reference_gene_annotation_source,
        params.reference_gene_annotation_assembly,
        params.reference_gene_annotation_version,
        params.reference_proteome_fasta_file,
        params.strategy,
        params.platform_tag,
        params.platform_unit_tag,
        params.library_tag,
        params.output_dir,
        params.exacto
    )
}
