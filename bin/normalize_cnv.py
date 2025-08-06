#!/usr/bin/env python3
"""
normalize_cnv.py
================
Aggregates exon coverage from all samples, performs z-score normalization
against reference cohort, and assigns copy-number states using fixed
thresholds.

Inputs
------
--coverage_dir   Directory containing per-sample ``*_coverage.tsv`` files.
--manifest       Manifest CSV (identical to previous modules).
--out            Path to output TSV summarising z-scores & CNV calls.

Output
------
TSV with columns:
    sample_id  gene  exon  zscore  copy_number
"""
import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

THRESHOLDS = [(-np.inf, -2.5, 0),
              (-2.5, -1.5, 1),
              (-1.5, 1.5, 2),
              (1.5, 2.5, 3),
              (2.5, np.inf, 4)]

def parse_args():
    p = argparse.ArgumentParser(description="Z-score normalisation of exon coverage and CN state inference.")
    p.add_argument("--coverage_dir", required=True, help="Directory with *_coverage.tsv files.")
    p.add_argument("--manifest", required=True, help="Manifest CSV to identify reference samples.")
    p.add_argument("--out", required=True, help="Output TSV with z-scores and copy numbers.")
    return p.parse_args()

def load_manifest(path: str) -> pd.DataFrame:
    return pd.read_csv(path)

def load_coverage(coverage_dir: str) -> pd.DataFrame:
    frames = []
    for cov_file in glob.glob(f"{coverage_dir.rstrip('/')}/*_coverage.tsv"):
        sample_id = Path(cov_file).stem.replace("_coverage", "")
        df = pd.read_csv(cov_file, sep="\t")
        df.insert(0, "sample_id", sample_id)
        frames.append(df)
    if not frames:
        sys.exit("No coverage files found. Check --coverage_dir argument.")
    return pd.concat(frames, ignore_index=True)

def compute_reference_stats(df: pd.DataFrame, reference_ids: list[str]):
    ref_df = df[df["sample_id"].isin(reference_ids)].copy()
    stats = ref_df.groupby(["gene", "exon"])["avg_depth"].agg(["mean", "std"]).reset_index()
    stats.rename(columns={"mean": "ref_mean", "std": "ref_std"}, inplace=True)
    return stats

def assign_copy_number(z: float) -> int:
    for low, high, cn in THRESHOLDS:
        if low < z <= high:
            return cn
    return 2  # fallback

def main():
    args = parse_args()
    manifest = load_manifest(args.manifest)
    coverage = load_coverage(args.coverage_dir)

    reference_ids = manifest.loc[manifest["type"].str.lower() == "reference", "sample_id"].tolist()
    if not reference_ids:
        sys.exit("No reference samples specified in manifest.")

    ref_stats = compute_reference_stats(coverage, reference_ids)
    merged = coverage.merge(ref_stats, on=["gene", "exon"], how="left")

    # Compute z-scores
    merged["zscore"] = (merged["avg_depth"] - merged["ref_mean"]) / merged["ref_std"]
    merged["copy_number"] = merged["zscore"].apply(assign_copy_number)

    out_cols = ["sample_id", "gene", "exon", "zscore", "copy_number"]
    merged[out_cols].to_csv(args.out, sep="\t", index=False)
    print(f"[normalize_cnv] Wrote {args.out}")

if __name__ == "__main__":
    main()