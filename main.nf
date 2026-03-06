#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    wes-pipeline: Whole Exome Sequencing Variant Calling Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github  : https://github.com/lorepaga1996/wes-pipeline
    Author  : lorepaga1996
    License : MIT
----------------------------------------------------------------------------------------
*/

log.info """\
    ╔══════════════════════════════════════════════╗
    ║           W E S - P I P E L I N E            ║
    ║   Whole Exome Sequencing Variant Calling     ║
    ╚══════════════════════════════════════════════╝

    samples     : ${params.samplesheet}
    genome      : ${params.genome}
    outdir      : ${params.outdir}
    """.stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC           } from './modules/qc/fastqc'
include { MULTIQC          } from './modules/qc/multiqc'
include { FASTP            } from './modules/trimming/fastp'
include { BWA_MEM2_INDEX   } from './modules/alignment/bwa_mem2_index'
include { BWA_MEM2_MEM     } from './modules/alignment/bwa_mem2_mem'
include { SAMTOOLS_SORT    } from './modules/alignment/samtools_sort'
include { SAMTOOLS_INDEX   } from './modules/alignment/samtools_index'
include { MARKDUPLICATES   } from './modules/alignment/markduplicates'
include { BASERECALIBRATOR } from './modules/variant_calling/baserecalibrator'
include { APPLYBQSR        } from './modules/variant_calling/applybqsr'
include { HAPLOTYPECALLER  } from './modules/variant_calling/haplotypecaller'
include { GENOTYPEGVCFS    } from './modules/variant_calling/genotypegvcfs'
include { VARIANTFILTRATION} from './modules/filtering/variantfiltration'
include { VEP              } from './modules/annotation/vep'
include { REPORT           } from './modules/filtering/report'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // ── 0. Input ────────────────────────────────────────────────────────────────
    ch_samplesheet = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [id: row.sample, sex: row.sex ?: 'unknown']
            [ meta, file(row.fastq_1), file(row.fastq_2) ]
        }

    ch_genome      = file(params.genome)
    ch_dbsnp       = file(params.dbsnp)
    ch_known_indels= file(params.known_indels)
    ch_vep_cache   = file(params.vep_cache)

    // ── 1. QC (pre-trim) ────────────────────────────────────────────────────────
    FASTQC(ch_samplesheet)

    // ── 2. Trimming ─────────────────────────────────────────────────────────────
    FASTP(ch_samplesheet)

    // ── 3. Alignment ────────────────────────────────────────────────────────────
    BWA_MEM2_INDEX(ch_genome)

    ch_trimmed = FASTP.out.reads
    BWA_MEM2_MEM(ch_trimmed, BWA_MEM2_INDEX.out.index)

    SAMTOOLS_SORT(BWA_MEM2_MEM.out.bam)
    MARKDUPLICATES(SAMTOOLS_SORT.out.bam)
    SAMTOOLS_INDEX(MARKDUPLICATES.out.bam)

    ch_bam_bai = MARKDUPLICATES.out.bam
        .join(SAMTOOLS_INDEX.out.bai)

    // ── 4. Base Quality Score Recalibration ─────────────────────────────────────
    BASERECALIBRATOR(ch_bam_bai, ch_genome, ch_dbsnp, ch_known_indels)
    APPLYBQSR(ch_bam_bai, BASERECALIBRATOR.out.table, ch_genome)

    ch_recal_bam = APPLYBQSR.out.bam
        .join(APPLYBQSR.out.bai)

    // ── 5. Variant Calling ──────────────────────────────────────────────────────
    HAPLOTYPECALLER(ch_recal_bam, ch_genome, file(params.intervals))
    GENOTYPEGVCFS(HAPLOTYPECALLER.out.gvcf, ch_genome)

    // ── 6. Variant Filtration ───────────────────────────────────────────────────
    VARIANTFILTRATION(GENOTYPEGVCFS.out.vcf, ch_genome)

    // ── 7. Annotation ───────────────────────────────────────────────────────────
    VEP(VARIANTFILTRATION.out.vcf, ch_vep_cache)

    // ── 8. MultiQC Report ───────────────────────────────────────────────────────
    ch_multiqc_files = FASTQC.out.zip
        .mix(FASTP.out.json)
        .mix(MARKDUPLICATES.out.metrics)
        .collect()

    MULTIQC(ch_multiqc_files)

    // ── 9. Final Report ─────────────────────────────────────────────────────────
    REPORT(VEP.out.vcf.collect(), MULTIQC.out.report)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    def status = workflow.success ? "✅ SUCCESS" : "❌ FAILED"
    log.info """
    Pipeline completed!
    Status   : ${status}
    Duration : ${workflow.duration}
    Results  : ${params.outdir}
    """
}
