# SMN1/SMN2 CNV Detection Pipeline (MVP)

This repository contains a minimal viable pipeline that detects copy-number variations (CNVs) in exons 7 and 8 of **SMN1** and **SMN2** using whole-exome sequencing (WES) data.

## Architecture

```
├── bin/                # Executable scripts
│   ├── run_pipeline.sh        # Master orchestrator
│   ├── extract_coverage.py    # Per-exon depth calculation
│   ├── allele_counts.py       # Discriminative SNP base counts
│   ├── normalize_cnv.py       # Z-score normalisation & CN calling
│   └── generate_report.py     # Human-readable reports
├── config/            # Static input files
│   ├── exons.bed              # GRCh38 exon coordinates
│   ├── smn_snps.tsv           # SMN1/2 discriminating SNPs
│   └── manifest.csv           # User-supplied sample sheet
└── results/           # Output (created on the fly)
```

## Quick Start

1. **Install dependencies**

   ```bash
   pip install -r requirements.txt
   ```

   * `samtools` ≥1.10 must be available on `$PATH` (or define `SAMTOOLS=/path/to/samtools`).

2. **Populate `config/manifest.csv`** with absolute BAM paths and label at least one sample as `reference`.

3. **Run the pipeline**

   ```bash
   bash bin/run_pipeline.sh
   ```

4. **Inspect results**

   * Per-sample coverage files: `results/coverage/*_coverage.tsv`
   * Allele counts: `results/allele_counts/*_alleles.tsv`
   * Z-scores & CN calls: `results/zscores.tsv`
   * Final human-readable reports: `results/reports/*_report.txt`

## Copy-Number State Thresholds

| Z-score range | Inferred CN | Interpretation |
| ------------- | ---------- | -------------- |
| ≤ −2.5        | 0          | Homozygous deletion |
| −2.5 – −1.5   | 1          | Heterozygous deletion |
| −1.5 – +1.5   | 2          | Normal diploid |
| +1.5 – +2.5   | 3          | Duplication |
| > +2.5        | 4+         | Multi-duplication |

## Notes

* Coordinates are approximate and should be adjusted against the latest genome build if high precision is required.
* The pipeline is intentionally lightweight and easy to extend. Feel free to replace shell calls to `samtools` with `pysam` or parallelise where appropriate.