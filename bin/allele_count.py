#!/usr/bin/env python3
import argparse, subprocess, csv

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', required=True)
    p.add_argument('--snps', required=True)
    p.add_argument('--out', required=True)
    return p.parse_args()

def read_manifest(manifest):
    samples = []
    with open(manifest) as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            sid, bam, stype = line.strip().split('\t')
            samples.append({'id': sid, 'bam': bam})
    return samples

def read_snps(snpfile):
    snps = []
    with open(snpfile) as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            chrom, pos, ref, alt, desc = line.strip().split('\t')
            snps.append({'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt})
    return snps

def count_alleles(bam, chrom, pos, ref, alt):
    cmd = [
        'samtools', 'mpileup', '-r', f'{chrom}:{pos}-{pos}', bam
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if not proc.stdout.strip(): return 0, 0
    fields = proc.stdout.strip().split('\t')
    if len(fields) < 5: return 0, 0
    bases = fields[4].upper()
    ref_count = bases.count(ref.upper())
    alt_count = bases.count(alt.upper())
    return ref_count, alt_count

def main():
    args = parse_args()
    samples = read_manifest(args.manifest)
    snps = read_snps(args.snps)
    with open(args.out, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['sample','chrom','pos','ref','alt','ref_count','alt_count'])
        for s in samples:
            for snp in snps:
                rc, ac = count_alleles(s['bam'], snp['chrom'], snp['pos'], snp['ref'], snp['alt'])
                writer.writerow([s['id'], snp['chrom'], snp['pos'], snp['ref'], snp['alt'], rc, ac])

if __name__ == '__main__':
    main()