#!/bin/bash
# Script to run ResFinder on all downloaded genomes
# Author: Vigneshwaran Muthuraman

# Activate ResFinder conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate resfinder

# Set directories
GENOME_DIR="data/genomes"
OUTPUT_DIR="results/resfinder"
DB_PATH="path/to/resfinder_db"  # resfinder database needs to be setup

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if database exists
if [ ! -d "$DB_PATH" ]; then
    echo "Error: ResFinder database not found at $DB_PATH"
    echo "Please download and set up the ResFinder database:"
    echo "  git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git $DB_PATH"
    echo "  cd $DB_PATH && python3 INSTALL.py"
    exit 1
fi

# Find all genome files
GENOMES=$(find "$GENOME_DIR" -name "*.fna" -o -name "*.fasta" -o -name "*.fa")
TOTAL=$(echo "$GENOMES" | wc -l)

echo "Running ResFinder on $TOTAL genomes..."

# Process each genome
COUNT=0
for GENOME in $GENOMES; do
    COUNT=$((COUNT + 1))
    
    # Extract genome ID from filename
    GENOME_ID=$(basename "$GENOME" | cut -d. -f1-2)
    echo "[$COUNT/$TOTAL] Processing $GENOME_ID..."
    
    # Create output directory for this genome
    mkdir -p "$OUTPUT_DIR/$GENOME_ID"
    
    # Run ResFinder
    echo "Running ResFinder on $GENOME_ID..."
    run_resfinder.py \
        -ifa "$GENOME" \
        -o "$OUTPUT_DIR/$GENOME_ID" \
        -s "acinetobacter" \
        -l 0.6 \
        -t 0.8 \
        --acquired \
        -db_res "$DB_PATH"
    
    # Check if analysis was successful
    if [ ! -f "$OUTPUT_DIR/$GENOME_ID/ResFinder_results_tab.txt" ]; then
        echo "Warning: ResFinder failed for $GENOME_ID"
    else
        echo "Successfully processed $GENOME_ID"
    fi
done

echo "ResFinder analysis complete."
