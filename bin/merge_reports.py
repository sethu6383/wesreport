#!/usr/bin/env python3
import sys
import os
import pandas as pd

manifest = sys.argv[1]
cnv_dir = sys.argv[2]
snp_dir = sys.argv[3]
out_dir = sys.argv[4]

os.makedirs(out_dir, exist_ok=True)

# Read manifest
df_manifest = pd.read_csv(manifest, sep='\t')
samples = df_manifest['sample_id'].tolist()

for sid in samples:
    cnv_file = os.path.join(cnv_dir, f"{sid}_cnv_report.tsv")
    snp_file = os.path.join(snp_dir, f"{sid}_snp_counts.tsv")
    if not os.path.exists(cnv_file):
        continue
    df_cnv = pd.read_csv(cnv_file, sep='\t')
    if os.path.exists(snp_file):
        df_snp = pd.read_csv(snp_file, sep='\t', names=['snp','chrom','pos','ref','alt','ref_count','alt_count'])
    else:
        df_snp = pd.DataFrame()
    with open(os.path.join(out_dir, f"{sid}_summary_report.tsv"), 'w') as out:
        out.write("# CNV Results\n")
        df_cnv.to_csv(out, sep='\t', index=False)
        out.write("\n# SNP Allele Counts\n")
        if not df_snp.empty:
            df_snp.to_csv(out, sep='\t', index=False)
        else:
            out.write("No SNP data available\n")