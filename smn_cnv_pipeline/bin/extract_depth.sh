#!/bin/bash

# extract_depth.sh - Extract read depth per exon from BAM files
# Usage: ./extract_depth.sh <sample_manifest> <bed_file> <output_dir>

set -euo pipefail

# Check arguments
if [ $# -ne 3 ]; then
    echo "Usage: $0 <sample_manifest> <bed_file> <output_dir>"
    echo "  sample_manifest: Path to sample manifest file"
    echo "  bed_file: Path to BED file with exon coordinates"
    echo "  output_dir: Output directory for depth files"
    exit 1
fi

MANIFEST="$1"
BED_FILE="$2"
OUTPUT_DIR="$3"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Log file
LOG_FILE="${OUTPUT_DIR}/depth_extraction.log"
echo "Starting depth extraction at $(date)" > "${LOG_FILE}"

# Process each sample
while IFS=$'\t' read -r sample_id bam_path sample_type population; do
    # Skip header and comment lines
    if [[ "$sample_id" =~ ^#.*$ ]] || [[ "$sample_id" == "sample_id" ]]; then
        continue
    fi
    
    echo "Processing sample: $sample_id" | tee -a "${LOG_FILE}"
    
    # Check if BAM file exists
    if [ ! -f "$bam_path" ]; then
        echo "WARNING: BAM file not found: $bam_path" | tee -a "${LOG_FILE}"
        continue
    fi
    
    # Check if BAM index exists
    if [ ! -f "${bam_path}.bai" ] && [ ! -f "${bam_path%.*}.bai" ]; then
        echo "WARNING: BAM index not found for: $bam_path" | tee -a "${LOG_FILE}"
        echo "Attempting to create index..." | tee -a "${LOG_FILE}"
        samtools index "$bam_path" 2>> "${LOG_FILE}" || {
            echo "ERROR: Failed to index BAM file: $bam_path" | tee -a "${LOG_FILE}"
            continue
        }
    fi
    
    # Output file for this sample
    depth_file="${OUTPUT_DIR}/${sample_id}_depth.txt"
    
    # Extract depth using samtools depth with BED file
    echo "Extracting depth for $sample_id..." | tee -a "${LOG_FILE}"
    samtools depth -b "$BED_FILE" -a "$bam_path" > "$depth_file" 2>> "${LOG_FILE}"
    
    # Check if depth extraction was successful
    if [ $? -eq 0 ] && [ -s "$depth_file" ]; then
        echo "Successfully extracted depth for $sample_id" | tee -a "${LOG_FILE}"
    else
        echo "ERROR: Depth extraction failed for $sample_id" | tee -a "${LOG_FILE}"
    fi
    
done < "$MANIFEST"

echo "Depth extraction completed at $(date)" | tee -a "${LOG_FILE}"
echo "Results saved in: $OUTPUT_DIR"