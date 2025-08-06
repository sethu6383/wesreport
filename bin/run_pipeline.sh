#!/usr/bin/env bash
set -euo pipefail

# Define sample information directly here
# Format: sample_id bam_path type
SAMPLES=(
  "REF1 /path/to/reference1.bam reference"
  "REF2 /path/to/reference2.bam reference"
  "TEST1 /path/to/test1.bam test"
)

BED="/workspace/config/smn_exons.bed"
SNPS="/workspace/config/smn_snps.tsv"
RESULTS="/workspace/results"

mkdir -p "$RESULTS"

# Write sample info to a temp file for downstream scripts
TMP_MANIFEST="$RESULTS/tmp_manifest.tsv"
echo -e "sample_id\tbam_path\ttype" > "$TMP_MANIFEST"
for entry in "${SAMPLES[@]}"; do
  echo -e "$entry" >> "$TMP_MANIFEST"
done

# 1. Extract per-exon coverage for all samples
python3 /workspace/bin/extract_coverage.py --manifest "$TMP_MANIFEST" --bed "$BED" --out "$RESULTS/coverage.tsv"

# 2. Allele-specific base counting at discriminating SNPs
python3 /workspace/bin/allele_count.py --manifest "$TMP_MANIFEST" --snps "$SNPS" --out "$RESULTS/allele_counts.tsv"

# 3. Normalize coverage and compute Z-scores
python3 /workspace/bin/normalize_coverage.py --manifest "$TMP_MANIFEST" --coverage "$RESULTS/coverage.tsv" --out "$RESULTS/zscores.tsv"

# 4. Call copy number states
python3 /workspace/bin/call_cn.py --zscores "$RESULTS/zscores.tsv" --out "$RESULTS/cn_calls.tsv"

# 5. Generate per-sample reports
python3 /workspace/bin/report.py --cn "$RESULTS/cn_calls.tsv" --allele "$RESULTS/allele_counts.tsv" --outdir "$RESULTS/reports"

echo "Pipeline complete. Results in $RESULTS."