#!/usr/bin/env python3
import argparse, csv
try:
    import numpy as np
except ImportError:
    print("Error: numpy is required. Install with 'pip install numpy'.")
    exit(1)

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', required=True)
    p.add_argument('--coverage', required=True)
    p.add_argument('--out', required=True)
    return p.parse_args()

def get_refs(manifest):
    refs = set()
    with open(manifest) as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            sid, _, stype = line.strip().split('\t')
            if stype == 'reference':
                refs.add(sid)
    return refs

def read_coverage(coverage):
    cov = {}
    with open(coverage) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            key = (row['sample'], row['exon'], row['gene'])
            cov[key] = float(row['mean_coverage'])
    return cov

def main():
    args = parse_args()
    refs = get_refs(args.manifest)
    cov = read_coverage(args.coverage)
    # Get all exons/genes
    exons = sorted(set((k[1], k[2]) for k in cov))
    samples = sorted(set(k[0] for k in cov))
    # Compute reference means and stds
    ref_stats = {}
    for exon, gene in exons:
        vals = [cov[(s, exon, gene)] for s in refs if (s, exon, gene) in cov]
        mu = np.mean(vals) if vals else 0
        sigma = np.std(vals, ddof=1) if len(vals) > 1 else 1
        ref_stats[(exon, gene)] = (mu, sigma)
    # Compute z-scores
    with open(args.out, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['sample','exon','gene','zscore'])
        for s in samples:
            for exon, gene in exons:
                v = cov.get((s, exon, gene), 0)
                mu, sigma = ref_stats[(exon, gene)]
                z = (v - mu) / sigma if sigma else 0
                writer.writerow([s, exon, gene, f'{z:.2f}'])

if __name__ == '__main__':
    main()