#!/bin/bash
# Script to run CRISPRCasFinder analysis on Acinetobacter genomes
# Author: Vigneshwaran Muthuraman

# Set up paths and directories
GENOME_DIR="data/genomes"
OUTPUT_DIR="results/crispr_analysis"
CONSOLIDATED_DIR="results/consolidated"
CRISPRCASFINDER_DIR="/path/to/CRISPRCasFinder"  # Update this to your actual path
CONDA_ENV_NAME="crisprcasfinder"  # Your conda environment name

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$CONSOLIDATED_DIR"

echo "=== CRISPRCasFinder Analysis Pipeline ==="
echo "Date: $(date)"

# Check if CRISPRCasFinder directory exists
if [ ! -d "$CRISPRCASFINDER_DIR" ]; then
    echo "Error: CRISPRCasFinder directory not found at $CRISPRCASFINDER_DIR"
    echo "Please update CRISPRCASFINDER_DIR variable in this script"
    echo "Or clone CRISPRCasFinder with: git clone https://github.com/dcouvin/CRISPRCasFinder.git"
    exit 1
fi

# Activate conda environment
echo "Activating conda environment: $CONDA_ENV_NAME"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$CONDA_ENV_NAME"

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate conda environment $CONDA_ENV_NAME"
    echo "Please create the environment first or update CONDA_ENV_NAME variable"
    exit 1
fi

# Change to CRISPRCasFinder directory (required for the tool to work)
echo "Changing to CRISPRCasFinder directory: $CRISPRCASFINDER_DIR"
cd "$CRISPRCASFINDER_DIR"

# Get absolute paths (since we're changing directories)
GENOME_DIR_ABS=$(realpath "../$GENOME_DIR")
OUTPUT_DIR_ABS=$(realpath "../$OUTPUT_DIR")

# Check if genome directory exists
if [ ! -d "$GENOME_DIR_ABS" ]; then
    echo "Error: Genome directory not found at $GENOME_DIR_ABS"
    exit 1
fi

# Find all genome files
GENOMES=$(find "$GENOME_DIR_ABS" -name "*.fna" -o -name "*.fasta" -o -name "*.fa")
GENOME_COUNT=$(echo "$GENOMES" | wc -l)

if [ $GENOME_COUNT -eq 0 ]; then
    echo "Error: No genome files found in $GENOME_DIR_ABS"
    exit 1
fi

echo "Found $GENOME_COUNT genome files"

# Combine all genomes into a single file
COMBINED_GENOMES="$OUTPUT_DIR_ABS/combined_genomes.fna"
echo "Combining genomes into: $COMBINED_GENOMES"

# Remove existing combined file if it exists
rm -f "$COMBINED_GENOMES"

# Combine genomes with proper headers
for GENOME in $GENOMES; do
    GENOME_ID=$(basename "$GENOME" | cut -d. -f1)
    echo "Adding $GENOME_ID to combined file..."
    
    # Add genome with modified header to include genome ID
    awk -v genome_id="$GENOME_ID" '
        /^>/ { 
            # Replace header with genome ID
            print ">" genome_id
            next
        }
        { print }
    ' "$GENOME" >> "$COMBINED_GENOMES"
done

echo "Combined genomes created: $COMBINED_GENOMES"
echo "Total sequences: $(grep -c "^>" "$COMBINED_GENOMES")"

# Run CRISPRCasFinder
echo "Running CRISPRCasFinder analysis..."
echo "Command: perl CRISPRCasFinder.pl -cas -keep -in $COMBINED_GENOMES -out $OUTPUT_DIR_ABS/crispr_results"

perl CRISPRCasFinder.pl \
    -cas \
    -keep \
    -in "$COMBINED_GENOMES" \
    -out "$OUTPUT_DIR_ABS/crispr_results"

# Check if analysis was successful
if [ $? -ne 0 ]; then
    echo "Error: CRISPRCasFinder analysis failed"
    exit 1
fi

echo "CRISPRCasFinder analysis completed successfully!"

# Process results
echo "Processing CRISPRCasFinder results..."

RESULTS_DIR="$OUTPUT_DIR_ABS/crispr_results"

# Check for output files
if [ -d "$RESULTS_DIR" ]; then
    echo "Results found in: $RESULTS_DIR"
    echo "Contents:"
    ls -la "$RESULTS_DIR"
    
    # Look for the summary file
    SUMMARY_FILE=""
    if [ -f "$RESULTS_DIR/result.json" ]; then
        SUMMARY_FILE="$RESULTS_DIR/result.json"
        echo "Found JSON result file: $SUMMARY_FILE"
    elif [ -f "$RESULTS_DIR/summary.txt" ]; then
        SUMMARY_FILE="$RESULTS_DIR/summary.txt"
        echo "Found summary file: $SUMMARY_FILE"
    fi
    
    # Look for TSV files
    TSV_FILES=$(find "$RESULTS_DIR" -name "*.tsv" 2>/dev/null)
    if [ -n "$TSV_FILES" ]; then
        echo "Found TSV files:"
        echo "$TSV_FILES"
    fi
    
    # Copy key result files to consolidated directory
    echo "Copying key results to consolidated directory..."
    
    # Copy any TSV files
    find "$RESULTS_DIR" -name "*.tsv" -exec cp {} "$CONSOLIDATED_DIR/" \;
    
    # Copy JSON if exists
    if [ -f "$RESULTS_DIR/result.json" ]; then
        cp "$RESULTS_DIR/result.json" "$CONSOLIDATED_DIR/crispr_cas_results.json"
    fi
    
    # Create a simple summary
    SUMMARY_OUTPUT="$CONSOLIDATED_DIR/crispr_cas_analysis_summary.txt"
    cat > "$SUMMARY_OUTPUT" << EOF
CRISPRCasFinder Analysis Summary
===============================
Date: $(date)
Input genomes: $GENOME_COUNT
Combined genome file: $COMBINED_GENOMES
Results directory: $RESULTS_DIR

Files created:
$(ls -1 "$CONSOLIDATED_DIR"/ | grep -E "(crispr|cas)" | sed 's/^/- /')

Next steps:
1. Check the results directory: $RESULTS_DIR
2. Look for TSV summary files
3. Run R analysis: Rscript code/2_analysis/8_crispr_cas_analysis.R
EOF
    
    echo "Summary created: $SUMMARY_OUTPUT"
    cat "$SUMMARY_OUTPUT"
    
else
    echo "Warning: Results directory not found at $RESULTS_DIR"
fi

# Return to original directory
cd - > /dev/null

echo ""
echo "=== CRISPRCasFinder Analysis Complete ==="
echo "Check results in: $OUTPUT_DIR"
echo "Consolidated files in: $CONSOLIDATED_DIR"
echo ""
echo "To run R analysis:"
echo "  Rscript code/2_analysis/8_crispr_cas_analysis.R"
