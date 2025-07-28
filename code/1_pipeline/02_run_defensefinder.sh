#!/bin/bash
# Script to run DefenseFinder on all downloaded genomes
# Author: Vigneshwaran Muthuraman

# Activate Defensefinder conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate defensefinder 

# Set directories
GENOME_DIR="data/genomes"
OUTPUT_DIR="results/defensefinder"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Find all genome files
GENOMES=$(find "$GENOME_DIR" -name "*.fna" -o -name "*.fasta" -o -name "*.fa")
TOTAL=$(echo "$GENOMES" | wc -l)

echo "Running DefenseFinder on $TOTAL genomes..."

# Process each genome
COUNT=0
for GENOME in $GENOMES; do
    COUNT=$((COUNT + 1))
    
    # Extract genome ID from filename
    GENOME_ID=$(basename "$GENOME" | cut -d. -f1-2) 
    echo "[$COUNT/$TOTAL] Processing $GENOME_ID..."
    
    # Create output directory for this genome
    mkdir -p "$OUTPUT_DIR/$GENOME_ID"
    
    # Run DefenseFinder
    # Note: This assumes DefenseFinder is installed and available in PATH
    # In a real environment, you might need to activate a conda environment first
    echo "Running defense-finder on $GENOME_ID..."
    defense-finder run -o "$OUTPUT_DIR/$GENOME_ID" "$GENOME"
    
    # Check if analysis was successful
    if [ ! -f "$OUTPUT_DIR/$GENOME_ID/${GENOME_ID}_defense_finder_systems.tsv" ]; then
        echo "Warning: DefenseFinder failed for $GENOME_ID"
    else
        echo "Successfully processed $GENOME_ID"
    fi
done

echo "DefenseFinder analysis complete."
