# SMN1/SMN2 CNV Detection MVP Pipeline

This pipeline detects copy number variations (CNVs) in the SMN1 and SMN2 genes (exons 7 and 8) from whole exome sequencing (WES) BAM files. It processes a cohort of samples (reference and test) and outputs per-sample CNV calls and supporting statistics.

## Folder Structure

- `bin/` — Executable scripts (bash, python)
- `config/` — Configuration files (BED, manifest, SNPs)
- `results/` — Output reports and intermediate files

## Pipeline Overview

1. **Input Preparation**
    - BED file: Exon coordinates for SMN1/2 (GRCh38)
    - Manifest: Sample BAM paths and type (reference/test)
    - SNPs: List of SMN1/2-discriminating SNPs
2. **Read Depth Extraction**
    - Uses `samtools depth` to extract per-exon coverage
    - Calculates average coverage per exon
3. **Allele-Specific Base Counting**
    - Uses pileup to count reads supporting each allele at discriminating SNPs
4. **Coverage Normalization & Z-score Calculation**
    - Normalizes coverage using reference samples
    - Computes per-exon Z-scores
5. **Copy Number Estimation**
    - Assigns CN states using Z-score thresholds
6. **Reporting**
    - Outputs per-sample reports with Z-scores, CN estimates, and allele counts

## Usage

1. Place your BED, manifest, and SNPs files in `config/`.
2. Place BAM files as referenced in the manifest.
3. Run the pipeline:

```bash
bash bin/run_pipeline.sh
```

Results will be in `results/`.

---

**Requirements:**
- bash, python3, samtools