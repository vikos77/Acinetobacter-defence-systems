#!/bin/bash
# Script to perform tblastn analysis of HMRG proteins against Acinetobacter genomes
# Author: Vigneshwaran Muthuraman

# Set directories
GENOME_DIR="data/genomes"
OUTPUT_DIR="results/hmrg_analysis"
HMRG_DIR="data/hmrg_database"
BLAST_RESULTS_DIR="$OUTPUT_DIR/blast_results"

# Create output directories
mkdir -p "$HMRG_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$BLAST_RESULTS_DIR"

# Download HMRG proteins from BacMet database
echo "Checking for BacMet HMRG database..."
if [ ! -f "$HMRG_DIR/BacMet2_EXP_database.fasta" ]; then
    echo "BacMet database file not found. You need to:"
    echo "1. Visit the BacMet website: https://bacmet.biomedicine.gu.se/download.html"
    echo "2. Download the BacMet2_EXP_database.fasta file"
    echo "3. Save it to $HMRG_DIR/BacMet2_EXP_database.fasta"
    echo "Once done, rerun this script."
    exit 1
fi

# Process the file - clean headers if needed
sed '/^>/ s/ /|/g' "$HMRG_DIR/BacMet2_EXP_database.fasta" > "$HMRG_DIR/hmrg_proteins.fasta"

# Combine all genomes into a single file for efficient BLAST
COMBINED_GENOME="$OUTPUT_DIR/all_genomes.fasta"
echo "Combining genomes for BLAST..."
cat $GENOME_DIR/*.f* > "$COMBINED_GENOME"

# Create the BLAST database using the genome data
echo "Formatting BLAST database..."
makeblastdb -in "$COMBINED_GENOME" -dbtype nucl -out results/databases/acinetobacter

# Run tblastn with stringent parameters for HMRG analysis
echo "Running tblastn against all genomes..."
tblastn -query "$HMRG_DIR/hmrg_proteins.fasta" \
       -db results/databases/acinetobacter \
       -out "$BLAST_RESULTS_DIR/HMRG_proteins_vs_acinetobacter.tblastn" \
       -evalue 0.005 \
       -qcov_hsp_perc 80 \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs qcovhsp" \
       -num_threads 8

echo "BLAST analysis completed. Results saved to $BLAST_RESULTS_DIR"
