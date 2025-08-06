#!/usr/bin/env bash
set -euo pipefail

MANIFEST=$1
BED=$2
OUTDIR=$3

mkdir -p "$OUTDIR"

while IFS=$'\t' read -r sample_id bam_path type; do
    [[ "$sample_id" == "sample_id" ]] && continue  # skip header
    echo "Processing $sample_id ($type)"
    while read -r chrom start end exon; do
        # samtools depth: -r region, -a for all positions
        region="$chrom:$((start+1))-$end"
        depth=$(samtools depth -r "$region" "$bam_path" | awk '{sum+=$3; n++} END {if(n>0) print sum/n; else print 0}')
        echo -e "$exon\t$depth" >> "$OUTDIR/${sample_id}_exon_coverage.tsv"
    done < "$BED"
done < "$MANIFEST"