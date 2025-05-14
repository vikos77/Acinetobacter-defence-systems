#!/bin/bash
# Script to run PADLOC on all downloaded genomes
# Author: Vigneshwaran Muthuraman

# Activate PADLOC conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate padloc


# Set directories
GENOME_DIR="data/genomes"
OUTPUT_DIR="results/padloc"

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
    BASE=$(basename "$GENOME") 
    GENOME_ID="${BASE%.*}" 
    echo "[$COUNT/$TOTAL] Processing $GENOME_ID..."
    
    # Create output directories for this genome
    mkdir -p "$OUTPUT_DIR/$GENOME_ID"
    
    # Run PADLOC
    echo "Running PADLOC for $GENOME_ID..."

        padloc --fna "$GENOME" --outdir "$OUTPUT_DIR/$GENOME_ID"
    
    # Check if analysis was successful
    if [ ! -f "$OUTPUT_DIR/$GENOME_ID/${GENOME_ID}_padloc.csv" ]; then
        echo "Warning: PADLOC failed for $GENOME_ID"
    else
        echo "Successfully processed $GENOME_ID"
    fi
done

echo "PADLOC analysis complete."
