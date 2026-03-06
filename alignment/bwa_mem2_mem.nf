process BWA_MEM2_INDEX {
    label 'process_high_memory'

    container 'quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5'

    publishDir "${params.outdir}/genome_index", mode: 'copy'

    input:
    path genome

    output:
    tuple path(genome), path("${genome}.*"), emit: index

    script:
    """
    bwa-mem2 index ${genome}
    """
}

process BWA_MEM2_MEM {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5'

    publishDir "${params.outdir}/alignment/raw", mode: 'copy'

    input:
    tuple val(meta), path(reads_1), path(reads_2)
    tuple path(genome), path(index)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam

    script:
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:${meta.id}"
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        ${genome} \\
        ${reads_1} ${reads_2} \\
    | samtools view -bS -@ ${task.cpus} -o ${meta.id}.bam
    """
}

process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    publishDir "${params.outdir}/alignment/sorted", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: bam

    script:
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -m ${(task.memory.toGiga() * 0.8 / task.cpus).intValue()}G \\
        -o ${meta.id}.sorted.bam \\
        ${bam}
    """
}

process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam}.bai"), emit: bai

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    """
}

process MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.outdir}/alignment/markdup", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.markdup.bam"), emit: bam
    path "${meta.id}.markdup.metrics",               emit: metrics

    script:
    """
    gatk MarkDuplicates \\
        --INPUT ${bam} \\
        --OUTPUT ${meta.id}.markdup.bam \\
        --METRICS_FILE ${meta.id}.markdup.metrics \\
        --REMOVE_DUPLICATES false \\
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
        --CREATE_INDEX false \\
        --java-options "-Xmx${(task.memory.toGiga() * 0.8).intValue()}g"
    """

    stub:
    """
    touch ${meta.id}.markdup.bam ${meta.id}.markdup.metrics
    """
}
