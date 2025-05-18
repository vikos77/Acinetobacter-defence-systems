#!/bin/bash
# Script to run AntiDefenseFinder on all downloaded genomes
# Author: Vigneshwaran Muthuraman

# Activate Defensefinder conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate defensefinder 

# Set directories
GENOME_DIR="data/genomes"
OUTPUT_DIR="results/antidefensefinder"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Find all genome files
GENOMES=$(find "$GENOME_DIR" -name "*.fna" -o -name "*.fasta" -o -name "*.fa")
TOTAL=$(echo "$GENOMES" | wc -l)

echo "Running AntiDefenseFinder on $TOTAL genomes..."

# Process each genome
COUNT=0
for GENOME in $GENOMES; do
    COUNT=$((COUNT + 1))
    
    # Extract genome ID from filename
    GENOME_ID=$(basename "$GENOME" | cut -d. -f1-2) 
    echo "[$COUNT/$TOTAL] Processing $GENOME_ID..."
    
    # Create output directory for this genome
    mkdir -p "$OUTPUT_DIR/$GENOME_ID"
    
    # Run DefenseFinder with AntiDefenseFinder only flag
    echo "Running anti-defense-finder on $GENOME_ID..."
    defense-finder run --antidefensefinder-only -o "$OUTPUT_DIR/$GENOME_ID" "$GENOME"
    
    # Check if analysis was successful
    if [ ! -f "$OUTPUT_DIR/$GENOME_ID/${GENOME_ID}_defense_finder_systems.tsv" ]; then
        echo "Warning: AntiDefenseFinder failed for $GENOME_ID"
    else
        echo "Successfully processed $GENOME_ID"
    fi
done

# Consolidate results
echo "Consolidating AntiDefenseFinder results..."
OUTPUT_CONSOLIDATED="results/consolidated/consolidated_antidefense_systems.tsv"
mkdir -p "results/consolidated"

# Create header line
echo -e "Genome_ID\tsys_id\ttype\tsubtype\tstart\tend\tstatus\tproduct\tgenes\tgenes_count\tcomments" > "$OUTPUT_CONSOLIDATED"

# Append data from all files
for FILE in $(find "$OUTPUT_DIR" -name "*_defense_finder_systems.tsv"); do
    GENOME_ID=$(basename "$(dirname "$FILE")")
    
    # Extract data, skip header, add genome ID
    awk -v genome="$GENOME_ID" 'NR>1 {print genome "\t" $0}' "$FILE" >> "$OUTPUT_CONSOLIDATED"
done

echo "AntiDefenseFinder analysis complete. Results consolidated to $OUTPUT_CONSOLIDATED"
