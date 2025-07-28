#!/bin/bash
# Script to download Acinetobacter genome sequences from NCBI
# Author: Vigneshwaran Muthuraman

# Set directories
OUTPUT_DIR="data/genomes"
ACCESSION_FILE="data/metadata/accession_list.txt"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if accession list exists
if [ ! -f "$ACCESSION_FILE" ]; then
    echo "Error: Accession list file not found at $ACCESSION_FILE"
    exit 1
fi

# Count total accessions
TOTAL=$(wc -l < "$ACCESSION_FILE")
echo "Downloading $TOTAL genomes..."

# Download each genome
COUNT=0
while read -r ACCESSION; do
    COUNT=$((COUNT + 1))
    echo "[$COUNT/$TOTAL] Downloading $ACCESSION..."
    
    # Use NCBI efetch to download genome
    ~/edirect/efetch -db nucleotide -id "$ACCESSION" -format fasta > "$OUTPUT_DIR/$ACCESSION.fna"
    
    # Check if download was successful
    if [ ! -s "$OUTPUT_DIR/$ACCESSION.fna" ]; then
        echo "Warning: Failed to download $ACCESSION"
    else
        echo "Successfully downloaded $ACCESSION"
    fi
    
    # Prevent overloading NCBI servers
    sleep 3
done < "$ACCESSION_FILE"

echo "Download complete. Downloaded $(ls -1 "$OUTPUT_DIR"/*.fna | wc -l) out of $TOTAL genomes."
