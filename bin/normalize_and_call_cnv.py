#!/usr/bin/env python3
import sys
import os
import pandas as pd
import numpy as np

manifest = sys.argv[1]
cov_dir = sys.argv[2]
out_dir = sys.argv[3]

os.makedirs(out_dir, exist_ok=True)

# Read manifest
df_manifest = pd.read_csv(manifest, sep='\t')
ref_samples = df_manifest[df_manifest['type'] == 'reference']['sample_id'].tolist()
test_samples = df_manifest[df_manifest['type'] == 'test']['sample_id'].tolist()

# Gather reference coverage
df_refs = []
for sid in ref_samples:
    f = os.path.join(cov_dir, f"{sid}_exon_coverage.tsv")
    df = pd.read_csv(f, sep='\t', names=['exon', 'depth'])
    df['sample_id'] = sid
    df_refs.append(df)
df_refs = pd.concat(df_refs)

# Compute reference mean and std per exon
ref_stats = df_refs.groupby('exon')['depth'].agg(['mean', 'std']).reset_index()

# Process test samples
for sid in test_samples:
    f = os.path.join(cov_dir, f"{sid}_exon_coverage.tsv")
    df = pd.read_csv(f, sep='\t', names=['exon', 'depth'])
    df = df.merge(ref_stats, on='exon')
    df['zscore'] = (df['depth'] - df['mean']) / df['std']
    # CNV calling
    def call_cnv(z):
        if z <= -2.5:
            return 0
        elif z <= -1.5:
            return 1
        elif z <= 1.5:
            return 2
        elif z <= 2.5:
            return 3
        else:
            return 4
    df['copy_number'] = df['zscore'].apply(call_cnv)
    df[['exon', 'depth', 'zscore', 'copy_number']].to_csv(os.path.join(out_dir, f"{sid}_cnv_report.tsv"), sep='\t', index=False)