process VEP {
    tag "$meta.id"
    label 'process_high'

    container 'ensemblorg/ensembl-vep:release_111'

    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  vep_cache

    output:
    tuple val(meta), path("${meta.id}.vep.vcf.gz"), emit: vcf
    path "${meta.id}.vep_summary.html",              emit: summary

    script:
    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${meta.id}.vep.vcf.gz \\
        --stats_file ${meta.id}.vep_summary.html \\
        --species ${params.vep_species} \\
        --assembly ${params.vep_assembly} \\
        --dir_cache ${vep_cache} \\
        --cache \\
        --offline \\
        --format vcf \\
        --vcf \\
        --compress_output bgzip \\
        --fork ${task.cpus} \\
        --everything \\
        --canonical \\
        --hgvs \\
        --symbol \\
        --numbers \\
        --domains \\
        --regulatory \\
        --af_gnomad \\
        --af_1kg \\
        --pubmed \\
        --variant_class \\
        --sift b \\
        --polyphen b \\
        --check_existing

    tabix -p vcf ${meta.id}.vep.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.vep.vcf.gz ${meta.id}.vep.vcf.gz.tbi ${meta.id}.vep_summary.html
    """
}
