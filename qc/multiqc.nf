process MULTIQC {
    label 'process_single'

    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'

    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    path(multiqc_files)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data

    script:
    def config = params.multiqc_config ? "--config ${params.multiqc_config}" : ''
    """
    multiqc \\
        --force \\
        ${config} \\
        --title "WES Pipeline QC Report" \\
        .
    """
}
