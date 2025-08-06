#!/usr/bin/env bash
set -euo pipefail

# Config
MANIFEST="/workspace/config/manifest.tsv"
BED="/workspace/config/smn_exons.bed"
SNPS="/workspace/config/smn_snps.tsv"
RESULTS="/workspace/results"

mkdir -p "$RESULTS"

# 1. Extract per-exon coverage for all samples
python3 /workspace/bin/extract_coverage.py --manifest "$MANIFEST" --bed "$BED" --out "$RESULTS/coverage.tsv"

# 2. Allele-specific base counting at discriminating SNPs
python3 /workspace/bin/allele_count.py --manifest "$MANIFEST" --snps "$SNPS" --out "$RESULTS/allele_counts.tsv"

# 3. Normalize coverage and compute Z-scores
python3 /workspace/bin/normalize_coverage.py --manifest "$MANIFEST" --coverage "$RESULTS/coverage.tsv" --out "$RESULTS/zscores.tsv"

# 4. Call copy number states
python3 /workspace/bin/call_cn.py --zscores "$RESULTS/zscores.tsv" --out "$RESULTS/cn_calls.tsv"

# 5. Generate per-sample reports
python3 /workspace/bin/report.py --cn "$RESULTS/cn_calls.tsv" --allele "$RESULTS/allele_counts.tsv" --outdir "$RESULTS/reports"

echo "Pipeline complete. Results in $RESULTS."