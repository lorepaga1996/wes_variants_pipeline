process BASERECALIBRATOR {
    tag "$meta.id"
    label 'process_medium'

    container 'broadinstitute/gatk:4.5.0.0'

    input:
    tuple val(meta), path(bam), path(bai)
    path  genome
    path  dbsnp
    path  known_indels

    output:
    tuple val(meta), path("${meta.id}.recal.table"), emit: table

    script:
    """
    gatk BaseRecalibrator \\
        -I ${bam} \\
        -R ${genome} \\
        --known-sites ${dbsnp} \\
        --known-sites ${known_indels} \\
        -O ${meta.id}.recal.table \\
        --java-options "-Xmx${(task.memory.toGiga() * 0.8).intValue()}g"
    """
}

process APPLYBQSR {
    tag "$meta.id"
    label 'process_medium'

    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.outdir}/alignment/recal", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(recal_table)
    path  genome

    output:
    tuple val(meta), path("${meta.id}.recal.bam"), emit: bam
    tuple val(meta), path("${meta.id}.recal.bai"), emit: bai

    script:
    """
    gatk ApplyBQSR \\
        -I ${bam} \\
        -R ${genome} \\
        --bqsr-recal-file ${recal_table} \\
        -O ${meta.id}.recal.bam \\
        --java-options "-Xmx${(task.memory.toGiga() * 0.8).intValue()}g"
    """
}

process HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high'

    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.outdir}/variant_calling/gvcf", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path  genome
    path  intervals

    output:
    tuple val(meta), path("${meta.id}.g.vcf.gz"), path("${meta.id}.g.vcf.gz.tbi"), emit: gvcf

    script:
    def intervals_opt = intervals ? "-L ${intervals}" : ''
    """
    gatk HaplotypeCaller \\
        -I ${bam} \\
        -R ${genome} \\
        -O ${meta.id}.g.vcf.gz \\
        -ERC GVCF \\
        ${intervals_opt} \\
        --native-pair-hmm-threads ${task.cpus} \\
        --java-options "-Xmx${(task.memory.toGiga() * 0.8).intValue()}g"
    """

    stub:
    """
    touch ${meta.id}.g.vcf.gz ${meta.id}.g.vcf.gz.tbi
    """
}

process GENOTYPEGVCFS {
    tag "$meta.id"
    label 'process_medium'

    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.outdir}/variant_calling/vcf", mode: 'copy'

    input:
    tuple val(meta), path(gvcf), path(tbi)
    path  genome

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), path("${meta.id}.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk GenotypeGVCFs \\
        -V ${gvcf} \\
        -R ${genome} \\
        -O ${meta.id}.vcf.gz \\
        --java-options "-Xmx${(task.memory.toGiga() * 0.8).intValue()}g"
    """
}
