#!/usr/bin/env nextflow

process runMoPepGenParseVEP {

    label 'mopepgen_parse'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tsv_file)
        path(reference_genome_fasta_file)
        path(reference_gtf_file)
        val(reference_source)
        val(variant_source)
        val(params_mopepgen_parsevep)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_${variant_source}_vep_mopepgen.gvf"), emit: f

    script:
        """
        moPepGen parseVEP \
            -i $tsv_file \
            --output-path ${sample_id}_${variant_source}_vep_mopepgen.gvf \
            --genome-fasta $reference_genome_fasta_file \
            --annotation-gtf $reference_gtf_file \
            --reference-source $reference_source \
            --source $variant_source \
            $params_mopepgen_parsevep

        if [ ! -f "${sample_id}_${variant_source}_vep_mopepgen.gvf" ]; then
            echo "##fileformat=VCFv4.2" > ${sample_id}_${variant_source}_vep_mopepgen.gvf
            echo "##parser=parseVEP" >> ${sample_id}_${variant_source}_vep_mopepgen.gvf
            echo "##reference_index=" >> ${sample_id}_${variant_source}_vep_mopepgen.gvf
            echo "##genome_fasta=" >> ${sample_id}_${variant_source}_vep_mopepgen.gvf
            echo "##annotation_gtf=" >> ${sample_id}_${variant_source}_vep_mopepgen.gvf
            echo "##source=${variant_source}" >> ${sample_id}_${variant_source}_vep_mopepgen.gvf
            echo "##CHROM=<Description='Gene ID'>" >> ${sample_id}_${variant_source}_vep_mopepgen.gvf
            echo '##INFO=<ID=TRANSCRIPT_ID,Number=1,Type=String,Description="Transcript ID">' >> ${sample_id}_${variant_source}_vep_mopepgen.gvf
            echo '##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">' >> ${sample_id}_${variant_source}_vep_mopepgen.gvf
            echo '##INFO=<ID=GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">' >> ${sample_id}_${variant_source}_vep_mopepgen.gvf
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> ${sample_id}_${variant_source}_vep_mopepgen.gvf
        fi
        """
}

process runMoPepGenParseREDItools2 {

    label 'mopepgen_parse'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tsv_file)
        path(reference_gtf_file)
        val(reference_source)
        val(params_mopepgen_parsereditools)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_reditools2_mopepgen.gvf"), emit: f

    script:
        """
        moPepGen parseREDItools \
            -i $tsv_file \
            --output-path ${sample_id}_reditools2_mopepgen.gvf \
            --annotation-gtf $reference_gtf_file \
            --reference-source $reference_source \
            --source reditools2 \
            $params_mopepgen_parsereditools

        if [ ! -f "${sample_id}_reditools2_mopepgen.gvf" ]; then
            echo "##fileformat=VCFv4.2" > ${sample_id}_reditools2_mopepgen.gvf
            echo "##parser=parseREDItools" >> ${sample_id}_reditools2_mopepgen.gvf
            echo "##reference_index=" >> ${sample_id}_reditools2_mopepgen.gvf
            echo "##genome_fasta=" >> ${sample_id}_reditools2_mopepgen.gvf
            echo "##annotation_gtf=" >> ${sample_id}_reditools2_mopepgen.gvf
            echo "##source=reditools2" >> ${sample_id}_reditools2_mopepgen.gvf
            echo "##CHROM=<Description='Gene ID'>" >> ${sample_id}_reditools2_mopepgen.gvf
            echo '##INFO=<ID=TRANSCRIPT_ID,Number=1,Type=String,Description="Transcript ID">' >> ${sample_id}_reditools2_mopepgen.gvf
            echo '##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">' >> ${sample_id}_reditools2_mopepgen.gvf
            echo '##INFO=<ID=GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">' >> ${sample_id}_reditools2_mopepgen.gvf
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> ${sample_id}_reditools2_mopepgen.gvf
        fi
        """
}

process runMoPepGenParseSTARFusion {

    label 'mopepgen_parse'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tsv_file)
        path(reference_genome_fasta_file)
        path(reference_genome_fasta_fai_file)
        path(reference_gtf_file)
        val(reference_source)
        val(params_mopepgen_parsestarfusion)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_starfusion_mopepgen.gvf"), emit: f

    script:
        """
        moPepGen parseSTARFusion \
            -i $tsv_file \
            --output-path ${sample_id}_starfusion_mopepgen.gvf \
            --genome-fasta $reference_genome_fasta_file \
            --annotation-gtf $reference_gtf_file \
            --reference-source $reference_source \
            --source starfusion \
            $params_mopepgen_parsestarfusion

        if [ ! -f "${sample_id}_starfusion_mopepgen.gvf" ]; then
            echo "##fileformat=VCFv4.2" > ${sample_id}_starfusion_mopepgen.gvf
            echo "##parser=parseSTARFusion" >> ${sample_id}_starfusion_mopepgen.gvf
            echo "##reference_index=" >> ${sample_id}_starfusion_mopepgen.gvf
            echo "##genome_fasta=" >> ${sample_id}_starfusion_mopepgen.gvf
            echo "##annotation_gtf=" >> ${sample_id}_starfusion_mopepgen.gvf
            echo "##source=starfusion" >> ${sample_id}_starfusion_mopepgen.gvf
            echo "##CHROM=<Description='Gene ID'>" >> ${sample_id}_starfusion_mopepgen.gvf
            echo '##INFO=<ID=TRANSCRIPT_ID,Number=1,Type=String,Description="Transcript ID">' >> ${sample_id}_starfusion_mopepgen.gvf
            echo '##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">' >> ${sample_id}_starfusion_mopepgen.gvf
            echo '##INFO=<ID=GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">' >> ${sample_id}_starfusion_mopepgen.gvf
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> ${sample_id}_starfusion_mopepgen.gvf
        fi
        """
}

process runMoPepGenParseArriba {

    label 'mopepgen_parse'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(tsv_file)
        path(reference_genome_fasta_file)
        path(reference_gtf_file)
        val(reference_source)
        val(params_mopepgen_parsearriba)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_arriba_mopepgen.gvf"), emit: f

    script:
        """
        moPepGen parseArriba \
            -i $tsv_file \
            --output-path ${sample_id}_arriba_mopepgen.gvf \
            --genome-fasta $reference_genome_fasta_file \
            --annotation-gtf $reference_gtf_file \
            --reference-source $reference_source \
            --source arriba \
            $params_mopepgen_parsearriba

        if [ ! -f "${sample_id}_arriba_mopepgen.gvf" ]; then
            echo "##fileformat=VCFv4.2" > ${sample_id}_arriba_mopepgen.gvf
            echo "##parser=parseArriba" >> ${sample_id}_arriba_mopepgen.gvf
            echo "##reference_index=" >> ${sample_id}_arriba_mopepgen.gvf
            echo "##genome_fasta=" >> ${sample_id}_arriba_mopepgen.gvf
            echo "##annotation_gtf=" >> ${sample_id}_arriba_mopepgen.gvf
            echo "##source=arriba" >> ${sample_id}_arriba_mopepgen.gvf
            echo "##CHROM=<Description='Gene ID'>" >> ${sample_id}_arriba_mopepgen.gvf
            echo '##INFO=<ID=TRANSCRIPT_ID,Number=1,Type=String,Description="Transcript ID">' >> ${sample_id}_arriba_mopepgen.gvf
            echo '##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">' >> ${sample_id}_arriba_mopepgen.gvf
            echo '##INFO=<ID=GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">' >> ${sample_id}_arriba_mopepgen.gvf
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> ${sample_id}_arriba_mopepgen.gvf
        fi
        """
}

process runMoPepGenParseRMATS {

    label 'mopepgen_parse'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(rmats_output_dir)
        path(reference_genome_fasta_file)
        path(reference_gtf_file)
        val(reference_source)
        val(params_mopepgen_parsermats)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_rmats_mopepgen.gvf"), emit: f

    script:
        """
        moPepGen parseRMATS \
            --se $rmats_output_dir/SE.MATS.JC.txt \
            --a5ss $rmats_output_dir/A5SS.MATS.JC.txt \
            --a3ss $rmats_output_dir/A3SS.MATS.JC.txt \
            --mxe $rmats_output_dir/MXE.MATS.JC.txt \
            --ri $rmats_output_dir/RI.MATS.JC.txt \
            --output-path ${sample_id}_rmats_mopepgen.gvf \
            --genome-fasta $reference_genome_fasta_file \
            --annotation-gtf $reference_gtf_file \
            --reference-source $reference_source \
            --source rmats \
            $params_mopepgen_parsermats

        if [ ! -f "${sample_id}_rmats_mopepgen.gvf" ]; then
            echo "##fileformat=VCFv4.2" > ${sample_id}_rmats_mopepgen.gvf
            echo "##parser=parseRMATS" >> ${sample_id}_rmats_mopepgen.gvf
            echo "##reference_index=" >> ${sample_id}_rmats_mopepgen.gvf
            echo "##genome_fasta=" >> ${sample_id}_rmats_mopepgen.gvf
            echo "##annotation_gtf=" >> ${sample_id}_rmats_mopepgen.gvf
            echo "##source=rmats" >> ${sample_id}_rmats_mopepgen.gvf
            echo "##CHROM=<Description='Gene ID'>" >> ${sample_id}_rmats_mopepgen.gvf
            echo '##INFO=<ID=TRANSCRIPT_ID,Number=1,Type=String,Description="Transcript ID">' >> ${sample_id}_rmats_mopepgen.gvf
            echo '##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">' >> ${sample_id}_rmats_mopepgen.gvf
            echo '##INFO=<ID=GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">' >> ${sample_id}_rmats_mopepgen.gvf
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> ${sample_id}_rmats_mopepgen.gvf
        fi
        """
}

process runMoPepGenParseCircExplorer2 {

    label 'mopepgen_parse'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(bed_file)
        path(reference_gtf_file)
        val(reference_source)
        val(params_mopepgen_parsecircexplorer)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_circexplorer2_mopepgen.gvf"), emit: f

    script:
        """
        cp $bed_file input.tsv

        moPepGen parseCIRCexplorer \
            -i input.tsv \
            --output-path ${sample_id}_circexplorer2_mopepgen.gvf \
            --annotation-gtf $reference_gtf_file \
            --reference-source $reference_source \
            --source circexplorer2 \
            $params_mopepgen_parsecircexplorer

        if [ ! -f "${sample_id}_circexplorer2_mopepgen.gvf" ]; then
            echo "##fileformat=VCFv4.2" > ${sample_id}_circexplorer2_mopepgen.gvf
            echo "##parser=parseCIRCexplorer" >> ${sample_id}_circexplorer2_mopepgen.gvf
            echo "##reference_index=" >> ${sample_id}_circexplorer2_mopepgen.gvf
            echo "##genome_fasta=" >> ${sample_id}_circexplorer2_mopepgen.gvf
            echo "##annotation_gtf=" >> ${sample_id}_circexplorer2_mopepgen.gvf
            echo "##source=circexplorer2" >> ${sample_id}_circexplorer2_mopepgen.gvf
            echo "##CHROM=<Description='Gene ID'>" >> ${sample_id}_circexplorer2_mopepgen.gvf
            echo '##INFO=<ID=TRANSCRIPT_ID,Number=1,Type=String,Description="Transcript ID">' >> ${sample_id}_circexplorer2_mopepgen.gvf
            echo '##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Gene Symbol">' >> ${sample_id}_circexplorer2_mopepgen.gvf
            echo '##INFO=<ID=GENOMIC_POSITION,Number=1,Type=String,Description="Genomic Position">' >> ${sample_id}_circexplorer2_mopepgen.gvf
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> ${sample_id}_circexplorer2_mopepgen.gvf
        fi
        """
}

process runMoPepGenCallVariant {

    label 'mopepgen_callvariant'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), val(gvf_files)
        path(reference_genome_fasta_file)
        path(reference_gtf_file)
        path(reference_proteome_fasta_file)
        val(params_mopepgen_callvariant)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_mopepgen_callvariant.fasta"), path("${sample_id}_mopepgen_callvariant_peptide_table.txt"), emit: f

    script:
        def input_files = gvf_files.findAll { it && it != 'null' && it != null }
        def input_files_str = input_files.join(' ')

        """
        moPepGen callVariant \
            -i ${input_files_str} \
            --output-path ${sample_id}_mopepgen_callvariant.fasta \
            --threads ${task.cpus} \
            --genome-fasta $reference_genome_fasta_file \
            --annotation-gtf $reference_gtf_file \
            --proteome-fasta $reference_proteome_fasta_file \
            $params_mopepgen_callvariant
        """
}
