#!/usr/bin/env bash
set -euo pipefail

# ---- USER CONFIGURABLE SECTION ----
# Define sample information here: sample_id bam_path type
SAMPLES=(
  "REF1 /path/to/reference1.bam reference"
  "REF2 /path/to/reference2.bam reference"
  "TEST1 /path/to/test1.bam test"
)
# Define exons (chrom start end exon gene)
EXONS=(
  "5 70247773 70247928 exon7 SMN1"
  "5 70249548 70249714 exon8 SMN1"
  "5 69372304 69372459 exon7 SMN2"
  "5 69374079 69374245 exon8 SMN2"
)
# Define SNPs (chrom pos ref alt description)
SNPS=(
  "5 70247773 C T SMN1/2_exon7_c.840C>T"
)
RESULTS="/workspace/results"
mkdir -p "$RESULTS"

# ---- COVERAGE EXTRACTION ----
echo -e "sample\texon\tgene\tmean_coverage" > "$RESULTS/coverage.tsv"
for entry in "${SAMPLES[@]}"; do
  read -r sid bam stype <<< "$entry"
  for exon in "${EXONS[@]}"; do
    read -r chrom start end exon_name gene <<< "$exon"
    mean_cov=$(samtools depth -r ${chrom}:${start}-${end} "$bam" | awk '{sum+=$3; n++} END {if(n>0) printf "%.2f", sum/n; else print 0}')
    echo -e "$sid\t$exon_name\t$gene\t$mean_cov" >> "$RESULTS/coverage.tsv"
  done
done

# ---- ALLELE COUNTING ----
echo -e "sample\tchrom\tpos\tref\talt\tref_count\talt_count" > "$RESULTS/allele_counts.tsv"
for entry in "${SAMPLES[@]}"; do
  read -r sid bam stype <<< "$entry"
  for snp in "${SNPS[@]}"; do
    read -r chrom pos ref alt desc <<< "$snp"
    pileup=$(samtools mpileup -r ${chrom}:${pos}-${pos} "$bam" 2>/dev/null)
    if [[ -z "$pileup" ]]; then refc=0; altc=0; else
      bases=$(echo "$pileup" | awk '{print toupper($5)}')
      refc=$(echo "$bases" | grep -o "$ref" | wc -l)
      altc=$(echo "$bases" | grep -o "$alt" | wc -l)
    fi
    echo -e "$sid\t$chrom\t$pos\t$ref\t$alt\t$refc\t$altc" >> "$RESULTS/allele_counts.tsv"
  done
done

# ---- NORMALIZATION & Z-SCORE ----
python3 - "$RESULTS/coverage.tsv" "$RESULTS/zscores.tsv" "${SAMPLES[@]}" <<'ENDPY'
import sys, csv
import numpy as np
covfile, out, *samples = sys.argv[1:]
refs = [s.split()[0] for s in samples if s.split()[2] == 'reference']
cov = {}
with open(covfile) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        key = (row['sample'], row['exon'], row['gene'])
        cov[key] = float(row['mean_coverage'])
exons = sorted(set((k[1], k[2]) for k in cov))
sample_ids = sorted(set(k[0] for k in cov))
ref_stats = {}
for exon, gene in exons:
    vals = [cov[(s, exon, gene)] for s in refs if (s, exon, gene) in cov]
    mu = np.mean(vals) if vals else 0
    sigma = np.std(vals, ddof=1) if len(vals) > 1 else 1
    ref_stats[(exon, gene)] = (mu, sigma)
with open(out, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['sample','exon','gene','zscore'])
    for s in sample_ids:
        for exon, gene in exons:
            v = cov.get((s, exon, gene), 0)
            mu, sigma = ref_stats[(exon, gene)]
            z = (v - mu) / sigma if sigma else 0
            writer.writerow([s, exon, gene, f'{z:.2f}'])
ENDPY

# ---- COPY NUMBER CALLING ----
echo -e "sample\texon\tgene\tzscore\tcn_call" > "$RESULTS/cn_calls.tsv"
awk 'NR==1{next} {z=$4+0; cn=(z<=-2.5)?0:(z<=-1.5)?1:(z<=1.5)?2:(z<=2.5)?3:4; print $1"\t"$2"\t"$3"\t"$4"\t"cn}' "$RESULTS/zscores.tsv" >> "$RESULTS/cn_calls.tsv"

# ---- REPORTING ----
REPORTDIR="$RESULTS/reports"
mkdir -p "$REPORTDIR"
for entry in "${SAMPLES[@]}"; do
  read -r sid bam stype <<< "$entry"
  rep="$REPORTDIR/${sid}_report.txt"
  echo "Report for sample: $sid" > "$rep"
  echo -e "\nCopy Number Calls (per exon):" >> "$rep"
  echo -e "exon\tgene\tzscore\tCN_call" >> "$rep"
  awk -v s="$sid" '$1==s {print $2"\t"$3"\t"$4"\t"$5}' "$RESULTS/cn_calls.tsv" >> "$rep"
  echo -e "\nAllele Counts (discriminating SNPs):" >> "$rep"
  echo -e "chrom\tpos\tref\talt\tref_count\talt_count" >> "$rep"
  awk -v s="$sid" '$1==s {print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' "$RESULTS/allele_counts.tsv" >> "$rep"
done

echo "Pipeline complete. Results in $RESULTS."