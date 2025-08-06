#!/usr/bin/env python3
"""
extract_coverage.py
===================
Module to compute average read depth over SMN1/SMN2 exons (7 & 8) for every
sample listed in a manifest. The script wraps `samtools depth` to obtain base
coverage and then summarises it per exon.

Inputs
------
--bed         BED file with columns: gene,exon,chrom,start,end (1-based, BED-like).
--manifest    CSV manifest with columns: sample_id,bam_path,type (type ∈ {reference,test}).
--out         Directory where per-sample coverage TSVs will be written.

Output
------
For each sample a TSV named ``<sample_id>_coverage.tsv`` with columns:
    gene   exon   chrom   start   end   avg_depth

This file is intended to be lightweight and easy to parse by subsequent
pipeline modules.
"""
import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path
from statistics import mean

import pandas as pd

SAMTOOLS_BIN = os.environ.get("SAMTOOLS", "samtools")

def parse_args():
    p = argparse.ArgumentParser(description="Compute average exon coverage using samtools depth.")
    p.add_argument("--bed", required=True, help="BED file with exon coordinates.")
    p.add_argument("--manifest", required=True, help="Manifest CSV with sample_id,bam_path,type.")
    p.add_argument("--out", required=True, help="Output directory for coverage TSVs.")
    p.add_argument("--threads", type=int, default=1, help="Unused placeholder for future parallelism.")
    return p.parse_args()

def load_bed(bed_path: str) -> pd.DataFrame:
    cols = ["gene", "exon", "chrom", "start", "end"]
    df = pd.read_csv(bed_path, sep="\t", names=cols, comment="#")
    return df

def load_manifest(manifest_path: str) -> pd.DataFrame:
    df = pd.read_csv(manifest_path)
    required_cols = {"sample_id", "bam_path", "type"}
    if not required_cols.issubset(df.columns):
        sys.exit(f"Manifest file must contain columns: {','.join(required_cols)}")
    return df

def samtools_depth_avg(bam: str, chrom: str, start: int, end: int) -> float:
    """Return average depth over region using samtools depth."""
    region = f"{chrom}:{start}-{end}"
    cmd = [SAMTOOLS_BIN, "depth", "-aa", "-r", region, bam]
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"samtools depth failed for {bam} region {region}: {e.stderr.strip()}")

    depths = []
    for line in proc.stdout.strip().splitlines():
        # format: chrom pos depth
        try:
            depth = int(line.split()[2])
        except IndexError:
            continue
        depths.append(depth)
    return mean(depths) if depths else 0.0

def process_sample(sample_row: pd.Series, bed_df: pd.DataFrame, out_dir: Path):
    sample_id = sample_row["sample_id"]
    bam = sample_row["bam_path"]
    out_file = out_dir / f"{sample_id}_coverage.tsv"

    records = []
    for _, row in bed_df.iterrows():
        avg_cov = samtools_depth_avg(bam, row.chrom, row.start, row.end)
        records.append({
            "gene": row.gene,
            "exon": row.exon,
            "chrom": row.chrom,
            "start": row.start,
            "end": row.end,
            "avg_depth": round(avg_cov, 2),
        })

    pd.DataFrame(records).to_csv(out_file, sep="\t", index=False)
    print(f"[extract_coverage] Wrote {out_file}")

def main():
    args = parse_args()
    bed_df = load_bed(args.bed)
    manifest_df = load_manifest(args.manifest)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    for _, sample_row in manifest_df.iterrows():
        process_sample(sample_row, bed_df, out_dir)

if __name__ == "__main__":
    main()