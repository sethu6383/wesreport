#!/usr/bin/env python3
"""
SMN Copy Number Variation (CNV) Analysis Pipeline

Improved version with better error handling, logging, configuration management,
and code organization for production use.

Author: Improved version
Date: 2024
"""

import os
import sys
import logging
import argparse
import json
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

import pysam
import pandas as pd
import numpy as np
from glob import glob


@dataclass
class CNVConfig:
    """Configuration class for CNV analysis parameters."""
    bam_folder: str
    ref_folder: str
    bed_file: str
    output_csv: str
    min_depth_threshold: float = 5.0
    min_ref_samples: int = 3
    z_score_thresholds: Dict[str, float] = None
    max_workers: int = 4
    
    def __post_init__(self):
        if self.z_score_thresholds is None:
            self.z_score_thresholds = {
                'homozygous_deletion': -2.5,
                'heterozygous_deletion': -1.2,
                'normal_low': -1.2,
                'normal_high': 1.2,
                'heterozygous_duplication': 2.5
            }


class CNVAnalyzer:
    """Main class for CNV analysis with improved error handling and logging."""
    
    def __init__(self, config: CNVConfig):
        self.config = config
        self.logger = self._setup_logging()
        self.ref_matrix = None
        self.ref_mean = None
        self.ref_std = None
        self.bed_df = None
        
    def _setup_logging(self) -> logging.Logger:
        """Setup logging configuration."""
        logger = logging.getLogger('CNVAnalyzer')
        logger.setLevel(logging.INFO)
        
        # Create console handler
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        
        # Create file handler
        log_dir = Path(self.config.output_csv).parent / "logs"
        log_dir.mkdir(exist_ok=True)
        file_handler = logging.FileHandler(log_dir / "cnv_analysis.log")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        
        return logger
    
    def validate_inputs(self) -> bool:
        """Validate all input files and directories exist."""
        errors = []
        
        # Check directories
        if not Path(self.config.bam_folder).exists():
            errors.append(f"BAM folder not found: {self.config.bam_folder}")
        if not Path(self.config.ref_folder).exists():
            errors.append(f"Reference folder not found: {self.config.ref_folder}")
            
        # Check BED file
        if not Path(self.config.bed_file).exists():
            errors.append(f"BED file not found: {self.config.bed_file}")
            
        # Check output directory is writable
        output_dir = Path(self.config.output_csv).parent
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except PermissionError:
            errors.append(f"Cannot create output directory: {output_dir}")
            
        if errors:
            for error in errors:
                self.logger.error(error)
            return False
            
        return True
    
    def load_and_validate_bams(self) -> Tuple[List[str], List[str]]:
        """Load and validate BAM files."""
        self.logger.info("Loading BAM files...")
        
        # Get BAM files
        ref_bams = sorted(glob(os.path.join(self.config.ref_folder, "*.bam")))
        sample_bams = sorted(glob(os.path.join(self.config.bam_folder, "*.bam")))
        
        # Remove reference BAMs from sample list
        ref_basenames = {os.path.basename(r) for r in ref_bams}
        sample_bams = [b for b in sample_bams if os.path.basename(b) not in ref_basenames]
        
        # Validate minimum reference samples
        if len(ref_bams) < self.config.min_ref_samples:
            raise ValueError(f"Insufficient reference samples: {len(ref_bams)} < {self.config.min_ref_samples}")
            
        # Validate BAM files
        invalid_bams = []
        for bam_list, name in [(ref_bams, "reference"), (sample_bams, "sample")]:
            for bam in bam_list:
                try:
                    with pysam.AlignmentFile(bam, "rb") as bam_file:
                        # Check if BAM is indexed
                        if not bam_file.has_index():
                            self.logger.warning(f"BAM file not indexed: {bam}")
                except Exception as e:
                    self.logger.error(f"Invalid {name} BAM file {bam}: {e}")
                    invalid_bams.append(bam)
        
        # Remove invalid BAMs
        ref_bams = [b for b in ref_bams if b not in invalid_bams]
        sample_bams = [b for b in sample_bams if b not in invalid_bams]
        
        self.logger.info(f"Loaded {len(ref_bams)} reference and {len(sample_bams)} sample BAMs")
        return ref_bams, sample_bams
    
    def load_bed_file(self) -> pd.DataFrame:
        """Load and validate BED file."""
        self.logger.info(f"Loading BED file: {self.config.bed_file}")
        
        try:
            bed_df = pd.read_csv(
                self.config.bed_file, 
                sep="\t", 
                header=None, 
                names=["chrom", "start", "end", "exon"],
                dtype={"chrom": str, "start": int, "end": int, "exon": str}
            )
            
            # Validate coordinates
            invalid_regions = bed_df[bed_df["start"] >= bed_df["end"]]
            if not invalid_regions.empty:
                self.logger.warning(f"Found {len(invalid_regions)} invalid regions with start >= end")
                bed_df = bed_df[bed_df["start"] < bed_df["end"]]
            
            # Create region strings
            bed_df["region"] = (bed_df["chrom"] + ":" + 
                              bed_df["start"].astype(str) + "-" + 
                              bed_df["end"].astype(str))
            
            # Add region length
            bed_df["length"] = bed_df["end"] - bed_df["start"]
            
            self.logger.info(f"Loaded {len(bed_df)} genomic regions")
            self.bed_df = bed_df
            return bed_df
            
        except Exception as e:
            self.logger.error(f"Error loading BED file: {e}")
            raise
    
    def calculate_depth_with_quality(self, bam_path: str, regions: pd.DataFrame, 
                                   min_mapq: int = 20) -> List[float]:
        """Calculate read depth with quality filtering."""
        try:
            bam = pysam.AlignmentFile(bam_path, "rb")
            depths = []
            
            for _, row in regions.iterrows():
                try:
                    # Count reads with quality filtering
                    count = 0
                    for read in bam.fetch(contig=row["chrom"], start=row["start"], end=row["end"]):
                        if (not read.is_unmapped and 
                            not read.is_duplicate and 
                            read.mapping_quality >= min_mapq):
                            count += 1
                    
                    depth = count / row["length"] if row["length"] > 0 else 0
                    depths.append(depth)
                    
                except Exception as e:
                    self.logger.warning(f"Error calculating depth for region {row['region']}: {e}")
                    depths.append(0.0)
            
            bam.close()
            return depths
            
        except Exception as e:
            self.logger.error(f"Error processing BAM file {bam_path}: {e}")
            raise
    
    def build_reference_matrix(self, ref_bams: List[str]) -> np.ndarray:
        """Build reference matrix with parallel processing."""
        self.logger.info("Building reference matrix...")
        
        ref_matrix = []
        
        # Use parallel processing for reference samples
        with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
            future_to_bam = {
                executor.submit(self.calculate_depth_with_quality, bam, self.bed_df): bam 
                for bam in ref_bams
            }
            
            for future in as_completed(future_to_bam):
                bam = future_to_bam[future]
                try:
                    depths = future.result()
                    ref_matrix.append(depths)
                    self.logger.info(f"Processed reference: {os.path.basename(bam)}")
                except Exception as e:
                    self.logger.error(f"Error processing reference {bam}: {e}")
        
        ref_matrix = np.array(ref_matrix)
        
        if ref_matrix.size == 0:
            raise ValueError("Reference matrix is empty")
        
        # Quality control: check for low-depth regions
        mean_depths = np.mean(ref_matrix, axis=0)
        low_depth_regions = np.where(mean_depths < self.config.min_depth_threshold)[0]
        
        if len(low_depth_regions) > 0:
            self.logger.warning(f"Found {len(low_depth_regions)} regions with low depth (< {self.config.min_depth_threshold})")
            for idx in low_depth_regions:
                region = self.bed_df.iloc[idx]["region"]
                self.logger.warning(f"Low depth region: {region} (depth: {mean_depths[idx]:.2f})")
        
        self.ref_matrix = ref_matrix
        return ref_matrix
    
    def calculate_normalization_parameters(self) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate mean and standard deviation for normalization."""
        self.logger.info("Calculating normalization parameters...")
        
        ref_mean = np.mean(self.ref_matrix, axis=0)
        ref_std = np.std(self.ref_matrix, axis=0, ddof=1)  # Use sample std
        
        # Handle zero standard deviation with small epsilon
        ref_std = np.where(ref_std == 0, 1e-6, ref_std)
        
        # Log statistics
        self.logger.info(f"Reference depth range: {ref_mean.min():.2f} - {ref_mean.max():.2f}")
        self.logger.info(f"Standard deviation range: {ref_std.min():.2f} - {ref_std.max():.2f}")
        
        self.ref_mean = ref_mean
        self.ref_std = ref_std
        
        return ref_mean, ref_std
    
    def classify_cnv(self, z_score: float) -> str:
        """Classify CNV based on Z-score with configurable thresholds."""
        thresholds = self.config.z_score_thresholds
        
        if z_score <= thresholds['homozygous_deletion']:
            return "Homozygous_Deletion"
        elif z_score <= thresholds['heterozygous_deletion']:
            return "Heterozygous_Deletion"
        elif thresholds['normal_low'] < z_score < thresholds['normal_high']:
            return "Normal"
        elif z_score < thresholds['heterozygous_duplication']:
            return "Heterozygous_Duplication"
        else:
            return "Homozygous_Duplication"
    
    def process_samples(self, sample_bams: List[str]) -> List[Dict]:
        """Process all sample BAMs and calculate CNV calls."""
        self.logger.info(f"Processing {len(sample_bams)} samples...")
        
        all_results = []
        
        for i, sample_bam in enumerate(sample_bams, 1):
            start_time = time.time()
            sample_name = os.path.basename(sample_bam).replace(".bam", "")
            
            try:
                # Calculate sample depths
                sample_depths = self.calculate_depth_with_quality(sample_bam, self.bed_df)
                
                # Calculate Z-scores
                z_scores = (np.array(sample_depths) - self.ref_mean) / self.ref_std
                
                # Classify CNVs
                classifications = [self.classify_cnv(z) for z in z_scores]
                
                # Store results
                for j, (_, exon) in enumerate(self.bed_df.iterrows()):
                    all_results.append({
                        "Sample": sample_name,
                        "Exon": exon["exon"],
                        "Region": exon["region"],
                        "Length": exon["length"],
                        "Raw_Depth": round(sample_depths[j], 2),
                        "Z_Score": round(z_scores[j], 2),
                        "CNV_Call": classifications[j],
                        "Reference_Mean": round(self.ref_mean[j], 2),
                        "Reference_Std": round(self.ref_std[j], 2)
                    })
                
                processing_time = time.time() - start_time
                self.logger.info(f"Processed sample {i}/{len(sample_bams)}: {sample_name} ({processing_time:.2f}s)")
                
            except Exception as e:
                self.logger.error(f"Error processing sample {sample_name}: {e}")
                continue
        
        return all_results
    
    def generate_summary_statistics(self, results_df: pd.DataFrame) -> Dict:
        """Generate summary statistics for the analysis."""
        summary = {
            "total_samples": results_df["Sample"].nunique(),
            "total_regions": results_df["Exon"].nunique(),
            "cnv_calls": results_df["CNV_Call"].value_counts().to_dict(),
            "samples_with_cnv": len(results_df[results_df["CNV_Call"] != "Normal"]["Sample"].unique()),
            "timestamp": pd.Timestamp.now().isoformat()
        }
        return summary
    
    def save_results(self, results: List[Dict]) -> None:
        """Save results with summary statistics."""
        self.logger.info(f"Saving results to: {self.config.output_csv}")
        
        # Create DataFrame
        results_df = pd.DataFrame(results)
        
        # Save main results
        results_df.to_csv(self.config.output_csv, index=False)
        
        # Generate and save summary
        summary = self.generate_summary_statistics(results_df)
        summary_file = Path(self.config.output_csv).with_suffix('.summary.json')
        
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        self.logger.info(f"Results saved: {len(results_df)} records")
        self.logger.info(f"Summary saved: {summary_file}")
        
        # Log CNV distribution
        cnv_counts = results_df["CNV_Call"].value_counts()
        self.logger.info("CNV distribution:")
        for cnv_type, count in cnv_counts.items():
            self.logger.info(f"  {cnv_type}: {count}")
    
    def run_analysis(self) -> None:
        """Main analysis pipeline."""
        start_time = time.time()
        self.logger.info("Starting CNV analysis pipeline...")
        
        try:
            # Validate inputs
            if not self.validate_inputs():
                raise ValueError("Input validation failed")
            
            # Load data
            ref_bams, sample_bams = self.load_and_validate_bams()
            self.load_bed_file()
            
            # Build reference
            self.build_reference_matrix(ref_bams)
            self.calculate_normalization_parameters()
            
            # Process samples
            results = self.process_samples(sample_bams)
            
            # Save results
            self.save_results(results)
            
            total_time = time.time() - start_time
            self.logger.info(f"✅ Analysis completed successfully in {total_time:.2f} seconds")
            
        except Exception as e:
            self.logger.error(f"Analysis failed: {e}")
            raise


def create_config_from_args() -> CNVConfig:
    """Create configuration from command line arguments."""
    parser = argparse.ArgumentParser(description="SMN CNV Analysis Pipeline")
    parser.add_argument("--bam-folder", required=True, help="Folder containing sample BAM files")
    parser.add_argument("--ref-folder", required=True, help="Folder containing reference BAM files")
    parser.add_argument("--bed-file", required=True, help="BED file with target regions")
    parser.add_argument("--output", required=True, help="Output CSV file")
    parser.add_argument("--min-depth", type=float, default=5.0, help="Minimum depth threshold")
    parser.add_argument("--min-refs", type=int, default=3, help="Minimum reference samples")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers")
    parser.add_argument("--config", help="JSON configuration file")
    
    args = parser.parse_args()
    
    if args.config and Path(args.config).exists():
        with open(args.config) as f:
            config_data = json.load(f)
        return CNVConfig(**config_data)
    else:
        return CNVConfig(
            bam_folder=args.bam_folder,
            ref_folder=args.ref_folder,
            bed_file=args.bed_file,
            output_csv=args.output,
            min_depth_threshold=args.min_depth,
            min_ref_samples=args.min_refs,
            max_workers=args.workers
        )


def main():
    """Main entry point."""
    try:
        config = create_config_from_args()
        analyzer = CNVAnalyzer(config)
        analyzer.run_analysis()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()