#!/bin/bash

# example_run.sh - Example of how to run the SMN CNV detection pipeline
# This is a demonstration script - modify paths for your actual data

echo "SMN CNV Detection Pipeline - Example Run"
echo "========================================"

# Set pipeline directory
PIPELINE_DIR="/workspace/smn_cnv_pipeline"
cd "$PIPELINE_DIR"

echo "1. First, edit the sample manifest file:"
echo "   nano config/sample_manifest.txt"
echo ""
echo "   Example content:"
echo "   sample_id	bam_path	sample_type	population"
echo "   REF001	/data/bams/ref001.bam	reference	EUR"
echo "   REF002	/data/bams/ref002.bam	reference	EUR"
echo "   TEST001	/data/bams/test001.bam	test	EUR"
echo ""

echo "2. Ensure your BAM files are indexed:"
echo "   samtools index /data/bams/ref001.bam"
echo "   samtools index /data/bams/ref002.bam"
echo "   samtools index /data/bams/test001.bam"
echo ""

echo "3. Run the pipeline:"
echo "   ./run_pipeline.sh"
echo ""
echo "   Or with custom directories:"
echo "   ./run_pipeline.sh --config ./config --results ./results"
echo ""

echo "4. Results will be available in:"
echo "   - results/reports/SAMPLE_ID/SAMPLE_ID_report.html (main reports)"
echo "   - results/cnv_calls/copy_numbers_gene_level.txt (summary)"
echo "   - results/pipeline_summary.txt (overall summary)"
echo ""

echo "5. For help:"
echo "   ./run_pipeline.sh --help"
echo ""

echo "Note: This is an example script. Modify the paths in config/sample_manifest.txt"
echo "to point to your actual BAM files before running the pipeline."