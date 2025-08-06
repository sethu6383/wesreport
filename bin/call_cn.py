#!/usr/bin/env python3
import argparse, csv

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--zscores', required=True)
    p.add_argument('--out', required=True)
    return p.parse_args()

def call_cn(z):
    z = float(z)
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

def main():
    args = parse_args()
    with open(args.zscores) as f, open(args.out, 'w', newline='') as out:
        reader = csv.DictReader(f, delimiter='\t')
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['sample','exon','gene','zscore','cn_call'])
        for row in reader:
            cn = call_cn(row['zscore'])
            writer.writerow([row['sample'], row['exon'], row['gene'], row['zscore'], cn])

if __name__ == '__main__':
    main()