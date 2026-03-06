# 🧬 wes-pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub last commit](https://img.shields.io/github/last-commit/<your-username>/wes-pipeline)](https://github.com/<your-username>/wes-pipeline)

A production-ready **Whole Exome Sequencing (WES) variant calling pipeline** built with [Nextflow DSL2](https://nextflow.io) and Docker. Implements GATK Best Practices for germline short variant discovery, from raw FASTQ to annotated VCF.

---

## Pipeline Overview

```
                     ┌─────────────────────────────────────────────────────┐
                     │               wes-pipeline workflow                 │
                     └─────────────────────────────────────────────────────┘

  FASTQ (R1 + R2)
       │
       ▼
  ┌─────────┐     ┌──────────────────────────────────────────────────────────┐
  │  FastQC │ ──► │                    QC REPORT                             │
  └─────────┘     │  FastQC (raw) ──► fastp ──► MultiQC                      │
       │          └──────────────────────────────────────────────────────────┘
       ▼
  ┌─────────┐
  │  fastp  │   Adapter trimming · Quality filtering · Length filtering
  └─────────┘
       │
       ▼
  ┌────────────┐
  │ BWA-MEM2   │   Read group tagging · Paired-end alignment (GRCh38)
  └────────────┘
       │
       ▼
  ┌──────────────────────┐
  │ Samtools sort + index│
  └──────────────────────┘
       │
       ▼
  ┌─────────────────┐
  │ MarkDuplicates  │   Optical + PCR duplicate flagging (GATK)
  └─────────────────┘
       │
       ▼
  ┌────────────────────────────────┐
  │   Base Quality Score Recalib.  │   BaseRecalibrator + ApplyBQSR
  └────────────────────────────────┘
       │
       ▼
  ┌──────────────────┐
  │ HaplotypeCaller  │   GVCF mode · per-sample variant calling
  └──────────────────┘
       │
       ▼
  ┌───────────────┐
  │ GenotypeGVCFs │   Joint genotyping
  └───────────────┘
       │
       ▼
  ┌────────────────────┐
  │ VariantFiltration  │   Hard filters: SNPs (QD, FS, MQ) · Indels (QD, FS)
  └────────────────────┘
       │
       ▼
  ┌─────────────────────────────────────────────────────────────┐
  │  VEP Annotation                                             │
  │  SIFT · PolyPhen · gnomAD AF · ClinVar · HGVS · Canonical   │
  └─────────────────────────────────────────────────────────────┘
       │
       ▼
  📊 HTML Summary Report + Annotated VCF
```

---

## Features

- **DSL2 modular design** — each tool is an independent, reusable module
- **Docker containers** — fully reproducible, no manual installations
- **GATK Best Practices** — BQSR, HaplotypeCaller GVCF mode, hard filtration
- **Comprehensive QC** — FastQC, fastp stats, MarkDuplicates metrics, MultiQC aggregation
- **VEP annotation** — population frequencies, functional predictions, clinical significance
- **Configurable profiles** — local, HPC (SLURM), cloud (AWS/Google)
- **Automatic retry** — handles transient failures with resource scaling

---

## Requirements

| Tool       | Version  |
|------------|----------|
| Nextflow   | ≥ 23.04  |
| Docker     | ≥ 20.10  |
| Java       | ≥ 11     |

---

## Quick Start

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
mv nextflow /usr/local/bin/
```

### 2. Clone the repository

```bash
git clone https://github.com/lorepaga1996/wes-pipeline.git
cd wes-pipeline
```

### 3. Prepare your samplesheet

Create a CSV file with the following columns:

```csv
sample,fastq_1,fastq_2,sex
PATIENT_01,/data/P01_R1.fastq.gz,/data/P01_R2.fastq.gz,female
PATIENT_02,/data/P02_R1.fastq.gz,/data/P02_R2.fastq.gz,male
```

### 4. Run the pipeline

```bash
nextflow run main.nf \
    --samplesheet samplesheet.csv \
    --genome /ref/GRCh38.fa \
    --intervals /ref/exome_targets.bed \
    --dbsnp /ref/dbsnp_146.hg38.vcf.gz \
    --known_indels /ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --vep_cache /ref/vep_cache \
    --outdir results \
    -profile docker
```

---

## Reference Files

Download GATK resource bundle (GRCh38):

```bash
# Genome
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta .

# dbSNP
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf .

# Known indels
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz .
```

Download VEP cache:

```bash
# GRCh38, Ensembl release 111
vep_install -a cf -s homo_sapiens -y GRCh38 --CACHEDIR /ref/vep_cache
```

---

## Configuration Profiles

| Profile  | Executor | Description                        |
|----------|----------|------------------------------------|
| `docker` | local    | Docker containers (default)        |
| `test`   | local    | Minimal test dataset               |
| `hpc`    | SLURM    | HPC cluster with SLURM scheduler   |
| `cloud`  | AWS/GCP  | Cloud execution                    |

---

## Output Structure

```
results/
├── qc/
│   ├── fastqc/             # Per-sample FastQC reports
│   └── multiqc/            # Aggregated MultiQC report
├── trimming/               # fastp trimmed FASTQs + stats
├── alignment/
│   ├── sorted/             # Coordinate-sorted BAMs
│   ├── markdup/            # Duplicate-marked BAMs + metrics
│   └── recal/              # BQSR-recalibrated BAMs
├── variant_calling/
│   ├── gvcf/               # Per-sample GVCFs
│   ├── vcf/                # Genotyped VCFs
│   └── filtered/           # Hard-filtered PASS VCFs
├── annotation/             # VEP-annotated VCFs + HTML summaries
└── report/                 # Final HTML summary report
```

---

## Variant Filtration Thresholds

### SNPs
| Filter          | Expression                 |
|-----------------|----------------------------|
| QD2             | `QD < 2.0`                 |
| FS60            | `FS > 60.0`                |
| MQ40            | `MQ < 40.0`                |
| MQRankSum-12.5  | `MQRankSum < -12.5`        |
| ReadPosRankSum-8| `ReadPosRankSum < -8.0`    |

### Indels
| Filter           | Expression                  |
|------------------|-----------------------------|
| QD2              | `QD < 2.0`                  |
| FS200            | `FS > 200.0`                |
| ReadPosRankSum-20| `ReadPosRankSum < -20.0`    |

---

## Tools & Containers

| Step                | Tool                   | Container                                    |
|---------------------|------------------------|----------------------------------------------|
| QC                  | FastQC 0.12.1          | `quay.io/biocontainers/fastqc:0.12.1`        |
| Trimming            | fastp 0.23.4           | `quay.io/biocontainers/fastp:0.23.4`         |
| Alignment           | BWA-MEM2 2.2.1         | `quay.io/biocontainers/bwa-mem2:2.2.1`       |
| Sort/Index          | SAMtools 1.19.2        | `quay.io/biocontainers/samtools:1.19.2`      |
| MarkDuplicates      | GATK 4.5.0             | `broadinstitute/gatk:4.5.0.0`                |
| BQSR / HC           | GATK 4.5.0             | `broadinstitute/gatk:4.5.0.0`                |
| Annotation          | VEP 111                | `ensemblorg/ensembl-vep:release_111`         |
| MultiQC             | MultiQC 1.21           | `quay.io/biocontainers/multiqc:1.21`         |

---

## References

1. Van der Auwera GA & O'Connor BD (2020). *Genomics in the Cloud*. O'Reilly Media.
2. DePristo MA et al. (2011). A framework for variation discovery and genotyping. *Nature Genetics*, 43, 491–498.
3. Di Tommaso P et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35, 316–319.
4. McLaren W et al. (2016). The Ensembl Variant Effect Predictor. *Genome Biology*, 17, 122.

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
