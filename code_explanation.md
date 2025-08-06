# SMN Copy Number Variation (CNV) Analysis Script

## Overview
This Python script performs **Copy Number Variation (CNV) analysis** for SMN (Survival Motor Neuron) genes using sequencing data in BAM format. It uses a Z-score based approach to detect deletions and duplications in SMN1 and SMN2 exons by comparing sample read depths to reference samples.

## What the Code Does

### 1. **Input Data Processing**
- **BAM Files**: Reads alignment files containing sequencing data
  - Reference samples (controls) from `/data/SMN/BAM/smn_z_score/refs/`
  - Test samples from `/data/SMN/BAM/`
- **BED File**: Defines genomic regions of interest (SMN1/SMN2 exons)
  - Format: chromosome, start position, end position, exon name

### 2. **Read Depth Calculation**
The `calculate_depth()` function:
- Opens each BAM file using `pysam`
- Counts aligned reads in each exon region
- Normalizes by region length to get reads per base pair
- Returns depth values for all exons

### 3. **Reference Baseline Establishment**
- Calculates read depths for all reference (control) samples
- Creates a reference matrix where:
  - Rows = reference samples
  - Columns = exons
- Computes mean and standard deviation for each exon across references

### 4. **Z-Score Normalization**
For each test sample:
```
Z-score = (Sample_Depth - Reference_Mean) / Reference_StdDev
```
This standardizes the data and makes it comparable across samples.

### 5. **CNV Classification**
Based on Z-scores, classifies each exon as:
- **Z ≤ -2.5**: Homozygous Deletion (both copies lost)
- **-2.5 < Z ≤ -1.2**: Heterozygous Deletion (one copy lost)
- **-1.2 < Z < 1.2**: Normal (typical copy number)
- **1.2 ≤ Z < 2.5**: Heterozygous Duplication (one extra copy)
- **Z ≥ 2.5**: Homozygous Duplication (multiple extra copies)

### 6. **Output Generation**
Creates a CSV file with results for each sample-exon combination:
- Sample name
- Exon identifier
- Genomic region
- Z-score
- Normalized read depth
- CNV classification

## Clinical Significance

**SMN genes** are crucial for:
- Motor neuron function
- Spinal Muscular Atrophy (SMA) diagnosis
- SMN1 deletions cause SMA
- SMN2 copy number affects disease severity

## Key Features

### Robust Statistics
- Handles zero standard deviation (prevents division by zero)
- Uses multiple reference samples for stability
- Normalizes for region length differences

### Quality Control
- Separates reference and test samples automatically
- Validates reference matrix isn't empty
- Creates output directories if needed

### Scalability
- Processes multiple samples in batch
- Handles variable numbers of exons
- Outputs structured data for downstream analysis

## Dependencies
- `pysam`: BAM file handling
- `pandas`: Data manipulation
- `numpy`: Numerical computations
- `glob`: File pattern matching

## Use Cases
1. **Clinical diagnostics**: SMA carrier screening
2. **Research**: Population studies of SMN copy numbers
3. **Quality control**: Validation of sequencing results
4. **Comparative analysis**: Before/after treatment studies