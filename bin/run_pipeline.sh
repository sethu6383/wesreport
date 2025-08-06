#!/usr/bin/env bash
set -euo pipefail

MANIFEST="/workspace/config/manifest.tsv"
BED="/workspace/config/smn_exons_grch38.bed"
COV_DIR="/workspace/results/coverage"
SNP_DIR="/workspace/results/snp_counts"
CNV_DIR="/workspace/results/cnv_calls"
REPORT_DIR="/workspace/results/reports"

mkdir -p "$COV_DIR" "$SNP_DIR" "$CNV_DIR" "$REPORT_DIR"

# 1. Extract exon coverage
echo "[1/4] Extracting exon coverage..."
/workspace/bin/extract_exon_coverage.sh "$MANIFEST" "$BED" "$COV_DIR"

# 2. Allele-specific base counting
echo "[2/4] Counting allele-specific bases..."
/workspace/bin/allele_specific_counts.sh "$MANIFEST" "$SNP_DIR"

# 3. Normalize and call CNV
echo "[3/4] Normalizing and calling CNV..."
/workspace/bin/normalize_and_call_cnv.py "$MANIFEST" "$COV_DIR" "$CNV_DIR"

# 4. Merge reports
echo "[4/4] Merging reports..."
/workspace/bin/merge_reports.py "$MANIFEST" "$CNV_DIR" "$SNP_DIR" "$REPORT_DIR"

echo "Pipeline complete. Reports in $REPORT_DIR"