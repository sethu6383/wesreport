#!/usr/bin/env python3

"""
normalize_coverage.py - Normalize coverage using reference samples and calculate Z-scores
Usage: python normalize_coverage.py <coverage_file> <sample_manifest> <output_file>
"""

import sys
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def read_sample_manifest(manifest_file):
    """Read sample manifest and return sample metadata."""
    samples = {}
    with open(manifest_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line or line.startswith('sample_id'):
                continue
            
            parts = line.split('\t')
            if len(parts) >= 3:
                sample_id, bam_path, sample_type = parts[:3]
                population = parts[3] if len(parts) > 3 else 'Unknown'
                samples[sample_id] = {
                    'bam_path': bam_path,
                    'sample_type': sample_type,
                    'population': population
                }
    return samples

def calculate_reference_stats(coverage_df, reference_samples):
    """Calculate mean and standard deviation for each exon using reference samples."""
    ref_stats = {}
    
    # Filter for reference samples only
    ref_df = coverage_df[coverage_df['sample_id'].isin(reference_samples)]
    
    if ref_df.empty:
        print("Warning: No reference samples found in coverage data!")
        return ref_stats
    
    # Calculate statistics for each exon
    exons = ref_df['exon'].unique()
    
    for exon in exons:
        exon_data = ref_df[ref_df['exon'] == exon]['avg_coverage']
        
        if len(exon_data) > 0:
            # Remove outliers using IQR method
            Q1 = exon_data.quantile(0.25)
            Q3 = exon_data.quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR
            
            # Filter outliers
            filtered_data = exon_data[(exon_data >= lower_bound) & (exon_data <= upper_bound)]
            
            if len(filtered_data) > 1:
                ref_stats[exon] = {
                    'mean': np.mean(filtered_data),
                    'std': np.std(filtered_data, ddof=1),
                    'median': np.median(filtered_data),
                    'n_samples': len(filtered_data),
                    'n_outliers': len(exon_data) - len(filtered_data),
                    'min': np.min(filtered_data),
                    'max': np.max(filtered_data)
                }
            else:
                # Fallback to all data if too few samples after filtering
                ref_stats[exon] = {
                    'mean': np.mean(exon_data),
                    'std': np.std(exon_data, ddof=1) if len(exon_data) > 1 else 1.0,
                    'median': np.median(exon_data),
                    'n_samples': len(exon_data),
                    'n_outliers': 0,
                    'min': np.min(exon_data),
                    'max': np.max(exon_data)
                }
        else:
            print(f"Warning: No coverage data for exon {exon} in reference samples")
    
    return ref_stats

def calculate_z_scores(coverage_df, ref_stats):
    """Calculate Z-scores for all samples using reference statistics."""
    z_score_data = []
    
    for _, row in coverage_df.iterrows():
        sample_id = row['sample_id']
        exon = row['exon']
        coverage = row['avg_coverage']
        
        if exon in ref_stats:
            ref_mean = ref_stats[exon]['mean']
            ref_std = ref_stats[exon]['std']
            
            # Avoid division by zero
            if ref_std > 0:
                z_score = (coverage - ref_mean) / ref_std
            else:
                z_score = 0.0
            
            z_score_data.append({
                'sample_id': sample_id,
                'exon': exon,
                'raw_coverage': coverage,
                'ref_mean': ref_mean,
                'ref_std': ref_std,
                'z_score': z_score,
                'ref_n_samples': ref_stats[exon]['n_samples']
            })
        else:
            print(f"Warning: No reference statistics for exon {exon}")
            z_score_data.append({
                'sample_id': sample_id,
                'exon': exon,
                'raw_coverage': coverage,
                'ref_mean': np.nan,
                'ref_std': np.nan,
                'z_score': np.nan,
                'ref_n_samples': 0
            })
    
    return pd.DataFrame(z_score_data)

def create_normalization_plots(z_scores_df, ref_stats, output_dir):
    """Create visualization plots for the normalization results."""
    plot_dir = Path(output_dir) / 'plots'
    plot_dir.mkdir(exist_ok=True)
    
    # Plot 1: Z-score distribution by exon
    plt.figure(figsize=(12, 8))
    exons = z_scores_df['exon'].unique()
    
    for i, exon in enumerate(sorted(exons)):
        plt.subplot(2, 2, i+1)
        exon_data = z_scores_df[z_scores_df['exon'] == exon]['z_score'].dropna()
        if len(exon_data) > 0:
            plt.hist(exon_data, bins=20, alpha=0.7, edgecolor='black')
            plt.axvline(x=0, color='red', linestyle='--', alpha=0.7)
            plt.title(f'{exon} Z-score Distribution')
            plt.xlabel('Z-score')
            plt.ylabel('Frequency')
    
    plt.tight_layout()
    plt.savefig(plot_dir / 'z_score_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Z-score heatmap
    pivot_z = z_scores_df.pivot(index='sample_id', columns='exon', values='z_score')
    
    plt.figure(figsize=(10, max(8, len(pivot_z) * 0.3)))
    sns.heatmap(pivot_z, cmap='RdBu_r', center=0, 
                cbar_kws={'label': 'Z-score'}, 
                fmt='.2f', linewidths=0.5)
    plt.title('Coverage Z-scores Heatmap')
    plt.ylabel('Sample ID')
    plt.xlabel('Exon')
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(plot_dir / 'z_score_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 3: Raw coverage vs Z-scores
    plt.figure(figsize=(15, 10))
    
    for i, exon in enumerate(sorted(exons)):
        plt.subplot(2, 2, i+1)
        exon_data = z_scores_df[z_scores_df['exon'] == exon]
        
        plt.scatter(exon_data['raw_coverage'], exon_data['z_score'], alpha=0.6)
        plt.axhline(y=0, color='red', linestyle='--', alpha=0.7)
        
        # Add reference mean line
        if exon in ref_stats:
            plt.axvline(x=ref_stats[exon]['mean'], color='green', 
                       linestyle='--', alpha=0.7, label='Ref Mean')
        
        plt.title(f'{exon}: Coverage vs Z-score')
        plt.xlabel('Raw Coverage')
        plt.ylabel('Z-score')
        plt.legend()
    
    plt.tight_layout()
    plt.savefig(plot_dir / 'coverage_vs_zscore.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved to: {plot_dir}")

def main():
    if len(sys.argv) != 4:
        print("Usage: python normalize_coverage.py <coverage_file> <sample_manifest> <output_file>")
        print("  coverage_file: Coverage data file (from calculate_coverage.py)")
        print("  sample_manifest: Sample manifest file")
        print("  output_file: Output file for normalized data")
        sys.exit(1)
    
    coverage_file = sys.argv[1]
    manifest_file = sys.argv[2]
    output_file = sys.argv[3]
    
    output_dir = Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Reading input files...")
    
    # Read coverage data
    try:
        coverage_df = pd.read_csv(coverage_file, sep='\t')
        print(f"Loaded coverage data for {len(coverage_df)} sample-exon combinations")
    except Exception as e:
        print(f"Error reading coverage file: {e}")
        sys.exit(1)
    
    # Read sample manifest
    samples_info = read_sample_manifest(manifest_file)
    
    # Identify reference samples
    reference_samples = [sid for sid, info in samples_info.items() 
                        if info['sample_type'] == 'reference']
    
    print(f"Found {len(reference_samples)} reference samples")
    
    if len(reference_samples) < 2:
        print("Warning: Very few reference samples. Results may be unreliable.")
    
    # Calculate reference statistics
    print("Calculating reference statistics...")
    ref_stats = calculate_reference_stats(coverage_df, reference_samples)
    
    # Calculate Z-scores
    print("Calculating Z-scores...")
    z_scores_df = calculate_z_scores(coverage_df, ref_stats)
    
    # Add sample type information
    z_scores_df['sample_type'] = z_scores_df['sample_id'].map(
        lambda x: samples_info.get(x, {}).get('sample_type', 'unknown')
    )
    z_scores_df['population'] = z_scores_df['sample_id'].map(
        lambda x: samples_info.get(x, {}).get('population', 'unknown')
    )
    
    # Save results
    z_scores_df.to_csv(output_file, index=False, sep='\t')
    
    # Save reference statistics
    ref_stats_file = output_file.replace('.txt', '_ref_stats.txt')
    ref_stats_df = pd.DataFrame.from_dict(ref_stats, orient='index').reset_index()
    ref_stats_df.rename(columns={'index': 'exon'}, inplace=True)
    ref_stats_df.to_csv(ref_stats_file, index=False, sep='\t')
    
    # Create plots
    try:
        create_normalization_plots(z_scores_df, ref_stats, output_dir)
    except Exception as e:
        print(f"Warning: Could not create plots: {e}")
    
    # Print summary
    print(f"\nNormalization completed!")
    print(f"Z-scores saved to: {output_file}")
    print(f"Reference statistics saved to: {ref_stats_file}")
    
    print(f"\nSummary statistics:")
    print(f"Total samples: {len(z_scores_df['sample_id'].unique())}")
    print(f"Reference samples: {len(reference_samples)}")
    print(f"Test samples: {len([s for s in samples_info.values() if s['sample_type'] == 'test'])}")
    
    print(f"\nZ-score summary by exon:")
    summary = z_scores_df.groupby('exon')['z_score'].agg(['count', 'mean', 'std', 'min', 'max']).round(3)
    print(summary)
    
    print(f"\nReference statistics:")
    print(ref_stats_df.round(3))

if __name__ == "__main__":
    main()