process VARIANTFILTRATION {
    tag "$meta.id"
    label 'process_medium'

    container 'broadinstitute/gatk:4.5.0.0'

    publishDir "${params.outdir}/variant_calling/filtered", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  genome

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi"), emit: vcf

    script:
    """
    # Split SNPs and indels, apply VQSR-like hard filters (suitable for single samples)

    # SNP filtration
    gatk SelectVariants -V ${vcf} -R ${genome} --select-type-to-include SNP -O snps.vcf.gz
    gatk VariantFiltration \\
        -V snps.vcf.gz \\
        -R ${genome} \\
        --filter-expression "QD < 2.0"    --filter-name "QD2" \\
        --filter-expression "FS > 60.0"   --filter-name "FS60" \\
        --filter-expression "MQ < 40.0"   --filter-name "MQ40" \\
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
        -O snps_filtered.vcf.gz

    # Indel filtration
    gatk SelectVariants -V ${vcf} -R ${genome} --select-type-to-include INDEL -O indels.vcf.gz
    gatk VariantFiltration \\
        -V indels.vcf.gz \\
        -R ${genome} \\
        --filter-expression "QD < 2.0"    --filter-name "QD2" \\
        --filter-expression "FS > 200.0"  --filter-name "FS200" \\
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \\
        -O indels_filtered.vcf.gz

    # Merge
    gatk MergeVcfs \\
        -I snps_filtered.vcf.gz \\
        -I indels_filtered.vcf.gz \\
        -O ${meta.id}.filtered.vcf.gz

    # Keep only PASS variants
    gatk SelectVariants \\
        -V ${meta.id}.filtered.vcf.gz \\
        -R ${genome} \\
        --exclude-filtered \\
        -O ${meta.id}.filtered.vcf.gz
    """
}

process REPORT {
    label 'process_single'

    container 'python:3.11-slim'

    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path vcfs
    path multiqc_report

    output:
    path "wes_summary_report.html", emit: html

    script:
    """
    #!/usr/bin/env python3
    import os, gzip, json
    from datetime import datetime

    # Count variants per sample
    summaries = []
    for vcf in "${vcfs}".split():
        sample = os.path.basename(vcf).replace('.filtered.vcf.gz','')
        count = 0
        snps = 0
        indels = 0
        try:
            with gzip.open(vcf, 'rt') as f:
                for line in f:
                    if not line.startswith('#'):
                        count += 1
                        ref, alt = line.split('\\t')[3], line.split('\\t')[4]
                        if len(ref) == len(alt) == 1:
                            snps += 1
                        else:
                            indels += 1
        except:
            pass
        summaries.append({'sample': sample, 'total': count, 'snps': snps, 'indels': indels})

    html = f'''<!DOCTYPE html>
<html>
<head><title>WES Pipeline Summary Report</title>
<style>
  body {{ font-family: monospace; background: #0d1117; color: #c9d1d9; padding: 40px; }}
  h1 {{ color: #58a6ff; }}
  table {{ border-collapse: collapse; width: 100%; margin-top: 20px; }}
  th, td {{ border: 1px solid #30363d; padding: 12px 16px; text-align: left; }}
  th {{ background: #161b22; color: #58a6ff; }}
  tr:hover {{ background: #161b22; }}
  .badge {{ background: #238636; color: white; padding: 2px 8px; border-radius: 12px; font-size: 12px; }}
</style>
</head>
<body>
  <h1>🧬 WES Pipeline Summary Report</h1>
  <p>Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}</p>
  <table>
    <tr><th>Sample</th><th>Total Variants</th><th>SNPs</th><th>Indels</th><th>Status</th></tr>
    {"".join(f"<tr><td>{s['sample']}</td><td>{s['total']}</td><td>{s['snps']}</td><td>{s['indels']}</td><td><span class='badge'>PASS</span></td></tr>" for s in summaries)}
  </table>
  <p style='margin-top:30px; color: #8b949e;'>Full QC report: <a href='../qc/multiqc/multiqc_report.html' style='color:#58a6ff;'>MultiQC</a></p>
</body>
</html>'''

    with open('wes_summary_report.html', 'w') as f:
        f.write(html)
    """
}
