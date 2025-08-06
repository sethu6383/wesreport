#!/usr/bin/env bash
set -euo pipefail

MANIFEST=$1
OUTDIR=$2

# Example SNPs: chrom, pos, ref, alt, name
SNP_LIST=(
    "chr5:70247773:C:T:SMN1_c.840C>T"
    "chr5:69372305:T:C:SMN2_c.840T>C"
)

mkdir -p "$OUTDIR"

for snp in "${SNP_LIST[@]}"; do
    IFS=":" read -r chrom pos ref alt name <<< "$snp"
    while IFS=$'\t' read -r sample_id bam_path type; do
        [[ "$sample_id" == "sample_id" ]] && continue
        counts=$(samtools mpileup -r "$chrom:$pos-$pos" -Q 20 -q 20 -f /path/to/reference.fa "$bam_path" \
            | awk -v ref="$ref" -v alt="$alt" -v pos="$pos" 'NR==1 {gsub("[.,]", "", $5); refc=gsub(ref, "", $5); altc=gsub(alt, "", $5); print refc, altc}')
        echo -e "$name\t$chrom\t$pos\t$ref\t$alt\t$counts" >> "$OUTDIR/${sample_id}_snp_counts.tsv"
    done < "$MANIFEST"
done