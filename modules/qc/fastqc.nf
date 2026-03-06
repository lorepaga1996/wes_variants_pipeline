process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads_1), path(reads_2)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${reads_1} ${reads_2}
    """

    stub:
    """
    touch ${meta.id}_R1_fastqc.html ${meta.id}_R1_fastqc.zip
    touch ${meta.id}_R2_fastqc.html ${meta.id}_R2_fastqc.zip
    """
}
