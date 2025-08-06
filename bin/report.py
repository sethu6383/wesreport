#!/usr/bin/env python3
import argparse, csv, os
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--cn', required=True)
    p.add_argument('--allele', required=True)
    p.add_argument('--outdir', required=True)
    return p.parse_args()

def read_tsv(path, keycol):
    d = defaultdict(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            d[row[keycol]].append(row)
    return d

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    cn = read_tsv(args.cn, 'sample')
    allele = read_tsv(args.allele, 'sample')
    for sample in set(cn) | set(allele):
        with open(os.path.join(args.outdir, f'{sample}_report.txt'), 'w') as out:
            out.write(f'Report for sample: {sample}\n')
            out.write('\nCopy Number Calls (per exon):\n')
            out.write('exon\tgene\tzscore\tCN_call\n')
            for row in cn.get(sample, []):
                out.write(f"{row['exon']}\t{row['gene']}\t{row['zscore']}\t{row['cn_call']}\n")
            out.write('\nAllele Counts (discriminating SNPs):\n')
            out.write('chrom\tpos\tref\talt\tref_count\talt_count\n')
            for row in allele.get(sample, []):
                out.write(f"{row['chrom']}\t{row['pos']}\t{row['ref']}\t{row['alt']}\t{row['ref_count']}\t{row['alt_count']}\n")
if __name__ == '__main__':
    main()