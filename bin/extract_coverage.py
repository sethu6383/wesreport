#!/usr/bin/env python3
import argparse, subprocess, csv, os

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', required=True)
    p.add_argument('--bed', required=True)
    p.add_argument('--out', required=True)
    return p.parse_args()

def read_bed(bedfile):
    exons = []
    with open(bedfile) as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            chrom, start, end, exon, gene = line.strip().split('\t')
            exons.append({'chrom': chrom, 'start': start, 'end': end, 'exon': exon, 'gene': gene})
    return exons

def read_manifest(manifest):
    samples = []
    with open(manifest) as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            sid, bam, stype = line.strip().split('\t')
            samples.append({'id': sid, 'bam': bam, 'type': stype})
    return samples

def mean_depth(bam, chrom, start, end):
    cmd = [
        'samtools', 'depth', '-r', f'{chrom}:{start}-{end}', bam
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    depths = [int(l.split()[2]) for l in proc.stdout.strip().split('\n') if l]
    return sum(depths)/len(depths) if depths else 0

def main():
    args = parse_args()
    exons = read_bed(args.bed)
    samples = read_manifest(args.manifest)
    with open(args.out, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['sample','exon','gene','mean_coverage'])
        for s in samples:
            for e in exons:
                cov = mean_depth(s['bam'], e['chrom'], e['start'], e['end'])
                writer.writerow([s['id'], e['exon'], e['gene'], f'{cov:.2f}'])

if __name__ == '__main__':
    main()