#!/usr/bin/env nextflow

process runNetMHCpan4 {

    label 'netmhcpan4'
    tag "${sample_id}"
    debug true

    publishDir(
        path: "${output_dir}/${sample_id}/",
        mode: 'copy'
    )

    input:
        tuple val(sample_id), path(fasta_file), val(mhc_alleles)
        val(netmhcpan_home_dir)
        val(netmhcpan)
        val(params_netmhcpan)
        val(output_dir)

    output:
        tuple val(sample_id), path("${sample_id}_netmhcpan4_outputs.txt"), emit: f

    script:
        """
        export NETMHCpan=$netmhcpan_home_dir
        $netmhcpan \
            -f $fasta_file \
            -a $mhc_alleles \
            -inptype 0 \
            $params_netmhcpan > ${sample_id}_netmhcpan4_outputs.txt
        """
}
