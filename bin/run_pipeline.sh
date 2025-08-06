#!/usr/bin/env bash
# run_pipeline.sh
# Master script that orchestrates SMN1/2 CNV detection pipeline.
# Usage: bin/run_pipeline.sh
set -euo pipefail

# Resolve script directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT="$( dirname "$DIR" )"

BED="$ROOT/config/exons.bed"
SNPS="$ROOT/config/smn_snps.tsv"
MANIFEST="$ROOT/config/manifest.csv"
RESULTS="$ROOT/results"
COV_DIR="$RESULTS/coverage"
ALLELE_DIR="$RESULTS/allele_counts"
REPORT_DIR="$RESULTS/reports"
ZSCORES_TSV="$RESULTS/zscores.tsv"

mkdir -p "$COV_DIR" "$ALLELE_DIR" "$REPORT_DIR"

echo "[pipeline] Extracting coverage..."
python "$DIR/extract_coverage.py" --bed "$BED" --manifest "$MANIFEST" --out "$COV_DIR"

echo "[pipeline] Counting alleles..."
python "$DIR/allele_counts.py" --snps "$SNPS" --manifest "$MANIFEST" --out "$ALLELE_DIR"

echo "[pipeline] Normalising coverage & calling CNVs..."
python "$DIR/normalize_cnv.py" --coverage_dir "$COV_DIR" --manifest "$MANIFEST" --out "$ZSCORES_TSV"

echo "[pipeline] Generating reports..."
python "$DIR/generate_report.py" --zscores "$ZSCORES_TSV" --allele_dir "$ALLELE_DIR" --out_dir "$REPORT_DIR"

echo "[pipeline] Completed. Reports are in $REPORT_DIR"