process FASTP {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'

    publishDir "${params.outdir}/trimming", mode: 'copy'

    input:
    tuple val(meta), path(reads_1), path(reads_2)

    output:
    tuple val(meta), path("${meta.id}_R1_trimmed.fastq.gz"), path("${meta.id}_R2_trimmed.fastq.gz"), emit: reads
    path "${meta.id}_fastp.json", emit: json
    path "${meta.id}_fastp.html", emit: html

    script:
    """
    fastp \\
        --in1 ${reads_1} \\
        --in2 ${reads_2} \\
        --out1 ${meta.id}_R1_trimmed.fastq.gz \\
        --out2 ${meta.id}_R2_trimmed.fastq.gz \\
        --json ${meta.id}_fastp.json \\
        --html ${meta.id}_fastp.html \\
        --detect_adapter_for_pe \\
        --correction \\
        --cut_right \\
        --cut_window_size 4 \\
        --cut_mean_quality 20 \\
        --length_required 50 \\
        --thread ${task.cpus}
    """

    stub:
    """
    touch ${meta.id}_R1_trimmed.fastq.gz ${meta.id}_R2_trimmed.fastq.gz
    touch ${meta.id}_fastp.json ${meta.id}_fastp.html
    """
}
