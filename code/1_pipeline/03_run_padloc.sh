#!/bin/bash
# Script to run PADLOC on all downloaded genomes using the wrapper script
# Author: Vigneshwaran Muthuraman

# Activate PADLOC conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate padloc

# Set directories
GENOME_DIR="data/genomes"
OUTPUT_DIR="results/padloc"
WRAPPER_SCRIPT="code/1_pipeline/padloc_wrapper.sh" 

# Make sure the wrapper script is executable
chmod +x "$WRAPPER_SCRIPT"

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR"

# Find all genome files
GENOMES=$(find "$GENOME_DIR" -name "*.fna" -o -name "*.fasta" -o -name "*.fa")
TOTAL=$(echo "$GENOMES" | wc -l)

echo "Running PADLOC on $TOTAL genomes..."

# Process each genome
COUNT=0
for GENOME in $GENOMES; do
    COUNT=$((COUNT + 1))
    
    # Extract genome ID from filename
    GENOME_ID=$(basename "$GENOME" | cut -d. -f1-2)
    echo "[$COUNT/$TOTAL] Processing $GENOME_ID..."
    
    # Create output directories for this genome
    mkdir -p "$OUTPUT_DIR/$GENOME_ID"
    
    # Run PADLOC using our wrapper script
    echo "Running PADLOC wrapper for $GENOME_ID..."
    
    # Run the wrapper script with output to the specific directory
    ./$WRAPPER_SCRIPT --cpu 8 --output "$OUTPUT_DIR/$GENOME_ID" "$GENOME"
    
    # Check if analysis was successful (note the changed filename pattern)
    if [ ! -f "$OUTPUT_DIR/$GENOME_ID/${GENOME_ID}_prodigal_padloc.csv" ]; then
        echo "Warning: PADLOC failed for $GENOME_ID"
    else
        echo "Successfully processed $GENOME_ID"
    fi
done

echo "PADLOC analysis complete."
