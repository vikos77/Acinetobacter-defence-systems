#!/bin/bash
# Script to perform tblastn analysis of IME proteins against Acinetobacter genomes
# Author: Vigneshwaran Muthuraman

# Set directories
GENOME_DIR="data/genomes"
OUTPUT_DIR="results/ime_analysis"
IME_DIR="data/ime_proteins/"
BLAST_RESULTS_DIR="$OUTPUT_DIR/blast_results"

# Create output directories
mkdir -p "$IME_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$BLAST_RESULTS_DIR"

# Download IME proteins from ICEberg (if available via direct download)
# Note: The ICEberg database might require manual download from the website
# In that case, you should place the downloaded files in $IME_DIR manually
# and comment out the wget commands

echo "Checking for IME protein database..."
if [ ! -f "$IME_DIR/IME_aa_all.fas.txt" ]; then
    echo "IME protein file not found. You need to:"
    echo "1. Visit the ICEberg website: https://db-mml.sjtu.edu.cn/ICEberg/"
    echo "2. Download the IME protein sequences"
    echo "3. Save them to $IME_DIR/ICEberg_proteins.faa"
    echo "Once done, rerun this script."
    exit 1
fi


# Combine all genomes into a single file for efficient BLAST
COMBINED_GENOME="$OUTPUT_DIR/all_genomes.fasta"
echo "Combining genomes for BLAST..."
cat $GENOME_DIR/*.f* > "$COMBINED_GENOME"

# Run tblastn
echo "Running tblastn against all genomes..."
tblastn -query "$IME_DIR/IME_aa_all.fas.txt" \
       -subject "$COMBINED_GENOME" \
       -out "$BLAST_RESULTS_DIR/IME_proteins_vs_acinetobacter.tblastn" \
       -evalue 1e-6 \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -num_threads 8

echo "BLAST analysis completed. Results saved to $BLAST_RESULTS_DIR"

# Optional: Parse BLAST results to create a summary file
echo "Creating summary of BLAST results..."
echo -e "Genome_ID\tIME_Type\tProtein_ID\tPercent_Identity\tE_value" > "$OUTPUT_DIR/ime_summary.tsv"

# Extract IME types from BLAST results
awk -F'\t' '$3 >= 80 && $11 <= 1e-6 {
    # Extract genome ID from subject sequence ID
    split($2, a, " ");
    genome_id = a[1];
    
    # Extract IME type from query sequence ID (assumes format like "ICEberg|123|protein_name")
    split($1, b, "|");
    if (length(b) >= 2) {
        ime_type = b[1] "|" b[2];
    } else {
        ime_type = $1;
    }
    
    # Extract protein ID
    protein_id = $1;
    
    # Output summary
    printf("%s\t%s\t%s\t%.2f\t%.2e\n", genome_id, ime_type, protein_id, $3, $11);
}' "$BLAST_RESULTS_DIR/IME_proteins_vs_acinetobacter.tblastn" >> "$OUTPUT_DIR/ime_summary.tsv"

echo "IME analysis complete. Summary saved to $OUTPUT_DIR/ime_summary.tsv"
