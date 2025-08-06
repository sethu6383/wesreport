#!/usr/bin/env python3
"""
allele_counts.py
================
Counts allele-specific read support at SMN1/SMN2 discriminating SNPs for each
sample in the manifest.

Inputs
------
--snps        TSV with columns: gene,chrom,pos,allele_sm1,allele_sm2
--manifest    CSV manifest with columns: sample_id,bam_path,type.
--out         Output directory for per-sample allele count TSVs.

Output
------
For each sample a TSV named ``<sample_id>_alleles.tsv`` with columns:
    gene  chrom  pos  allele  count
Where *allele* is one of the two specified discriminating alleles.
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

SAMTOOLS_BIN = os.environ.get("SAMTOOLS", "samtools")

def parse_args():
    p = argparse.ArgumentParser(description="Count allele support at discriminating SNPs using samtools mpileup.")
    p.add_argument("--snps", required=True, help="TSV with discriminating SNPs.")
    p.add_argument("--manifest", required=True, help="Manifest CSV with sample_id,bam_path,type.")
    p.add_argument("--out", required=True, help="Output directory for allele counts.")
    return p.parse_args()

def load_snps(snp_path: str) -> pd.DataFrame:
    cols = ["gene", "chrom", "pos", "allele_sm1", "allele_sm2"]
    df = pd.read_csv(snp_path, sep="\t", names=cols, comment="#")
    return df

def load_manifest(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    required_cols = {"sample_id", "bam_path"}
    if not required_cols.issubset(df.columns):
        sys.exit(f"Manifest missing columns: {','.join(required_cols)}")
    return df

def mpileup_counts(bam: str, chrom: str, pos: int) -> str:
    """Return pileup string for a single genomic position."""
    region = f"{chrom}:{pos}-{pos}"
    cmd = [SAMTOOLS_BIN, "mpileup", "-aa", "-Q", "0", "-r", region, bam]
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"samtools mpileup failed: {e.stderr.strip()}")
    return proc.stdout.strip()

def count_bases(pileup_line: str, allele1: str, allele2: str) -> tuple[int, int]:
    """Count how many bases match allele1 and allele2 in the pileup string."""
    if not pileup_line:
        return 0, 0
    fields = pileup_line.split("\t")
    if len(fields) < 5:
        return 0, 0
    bases = fields[4]
    # Pileup encoding is complex; lower-/uppercase for forward/reverse read, '.' ',' for matches etc.
    # We'll simplify by converting to upper case and counting letters.
    bases_clean = bases.replace(".", fields[2].upper()).replace(",", fields[2].upper()).upper()
    count1 = bases_clean.count(allele1.upper())
    count2 = bases_clean.count(allele2.upper())
    return count1, count2

def process_sample(sample_row: pd.Series, snp_df: pd.DataFrame, out_dir: Path):
    sample_id = sample_row["sample_id"]
    bam = sample_row["bam_path"]

    records = []
    for _, snp in snp_df.iterrows():
        pile = mpileup_counts(bam, snp.chrom, snp.pos)
        c1, c2 = count_bases(pile, snp.allele_sm1, snp.allele_sm2)
        records.append({
            "gene": snp.gene,
            "chrom": snp.chrom,
            "pos": snp.pos,
            "allele": snp.allele_sm1,
            "count": c1,
        })
        records.append({
            "gene": snp.gene,
            "chrom": snp.chrom,
            "pos": snp.pos,
            "allele": snp.allele_sm2,
            "count": c2,
        })

    out_file = out_dir / f"{sample_id}_alleles.tsv"
    pd.DataFrame(records).to_csv(out_file, sep="\t", index=False)
    print(f"[allele_counts] Wrote {out_file}")

def main():
    args = parse_args()
    snp_df = load_snps(args.snps)
    manifest_df = load_manifest(args.manifest)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    for _, sample_row in manifest_df.iterrows():
        process_sample(sample_row, snp_df, out_dir)

if __name__ == "__main__":
    main()