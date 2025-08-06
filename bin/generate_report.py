#!/usr/bin/env python3
"""
generate_report.py
==================
Combines z-score/CN calls with allele counts to produce per-sample summary
reports.

Inputs
------
--zscores     Output TSV from `normalize_cnv.py`.
--allele_dir  Directory with ``*_alleles.tsv`` files from `allele_counts.py`.
--out_dir     Directory where the per-sample reports will be saved (plain text).
"""
import argparse
import glob
from pathlib import Path

import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Generate per-sample CNV reports for SMN1/SMN2.")
    p.add_argument("--zscores", required=True, help="TSV with z-scores & CN calls.")
    p.add_argument("--allele_dir", required=True, help="Directory containing *_alleles.tsv files.")
    p.add_argument("--out_dir", required=True, help="Output directory for reports.")
    return p.parse_args()

def load_allele_counts(allele_dir: str) -> dict[str, pd.DataFrame]:
    mapping = {}
    for f in glob.glob(f"{allele_dir.rstrip('/') }/*_alleles.tsv"):
        sample_id = Path(f).stem.replace("_alleles", "")
        mapping[sample_id] = pd.read_csv(f, sep="\t")
    return mapping

def write_report(sample_id: str, z_df: pd.DataFrame, allele_df: pd.DataFrame, out_dir: Path):
    lines = [f"Sample: {sample_id}", "=" * 40]
    for (_, row) in z_df.iterrows():
        gene_exon = f"{row['gene']} exon{row['exon']}"
        lines.append(f"{gene_exon:15s}  Z-score: {row['zscore']:.2f}  CN: {row['copy_number']}")
    lines.append("\nAllele counts (discriminating SNPs):")
    if allele_df is not None:
        for _, r in allele_df.iterrows():
            lines.append(f"{r['gene']} {r['chrom']}:{r['pos']} {r['allele']} => {r['count']}")
    else:
        lines.append("No allele information available.")
    report_path = out_dir / f"{sample_id}_report.txt"
    report_path.write_text("\n".join(lines))
    print(f"[generate_report] Wrote {report_path}")

def main():
    args = parse_args()
    zscores_df = pd.read_csv(args.zscores, sep="\t")
    allele_map = load_allele_counts(args.allele_dir)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for sample_id, group in zscores_df.groupby("sample_id"):
        allele_df = allele_map.get(sample_id)
        write_report(sample_id, group, allele_df, out_dir)

if __name__ == "__main__":
    main()