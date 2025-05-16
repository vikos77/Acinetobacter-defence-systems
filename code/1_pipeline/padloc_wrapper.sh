#!/bin/bash
# padloc_wrapper.sh - A comprehensive wrapper for PADLOC analysis

# Display help information
show_help() {
    echo "Usage: $0 [options] genome.fna"
    echo "A wrapper script that properly formats Prodigal output for PADLOC analysis."
    echo ""
    echo "Options:"
    echo "  -c, --cpu NUM       Number of CPU cores to use (default: 1)"
    echo "  -f, --force         Overwrite existing files"
    echo "  -o, --output DIR    Output directory (default: current directory)"
    echo "  -h, --help          Display this help message"
    echo ""
    echo "Example: $0 --cpu 8 genome.fna"
    echo ""
    exit 1
}

# Default values
CPU=1
FORCE=0
OUTPUT_DIR="."

# Parse command line arguments
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -c|--cpu)
            CPU="$2"
            shift 2
            ;;
        -f|--force)
            FORCE=1
            shift
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

# Check if a genome file was provided
if [ $# -ne 1 ]; then
    echo "Error: You must specify a genome file."
    show_help
fi

GENOME="$1"
BASENAME=$(basename "$GENOME" .fna)

# Make sure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Use temp directory for intermediate files
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT

# Step 1: Run Prodigal to generate protein sequences and GFF
echo "[$(date +%T)] Running Prodigal for $GENOME..."
if [ $FORCE -eq 1 ] || [ ! -f "${TEMP_DIR}/${BASENAME}_prodigal.faa" ] || [ ! -f "${TEMP_DIR}/${BASENAME}_prodigal.gff" ]; then
    prodigal -i "$GENOME" -a "${TEMP_DIR}/${BASENAME}_prodigal.faa" -f gff -o "${TEMP_DIR}/${BASENAME}_prodigal.gff"
    if [ $? -ne 0 ]; then
        echo "Error: Prodigal failed. Exiting."
        exit 1
    fi
else
    echo "[$(date +%T)] Prodigal output files already exist. Use --force to overwrite."
fi

# Step 2: Extract protein IDs from the FASTA file
echo "[$(date +%T)] Extracting protein IDs from FASTA..."
grep ">" "${TEMP_DIR}/${BASENAME}_prodigal.faa" | sed 's/>//' > "${TEMP_DIR}/protein_ids.txt"

# Step 3: Create a properly formatted GFF file
echo "[$(date +%T)] Creating properly formatted GFF file..."
awk -v OFS="\t" '
BEGIN {
    # Read protein IDs from file
    while ((getline id < "'${TEMP_DIR}'/protein_ids.txt") > 0) {
        protein_ids[++count] = id;
    }
    close("'${TEMP_DIR}'/protein_ids.txt");
    current_id = 1;
}
/^#/ {
    print $0;  # Print comment lines unchanged
    next;
}
$3 == "CDS" {
    # Add ID attribute that exactly matches the FASTA header
    $9 = "ID=" protein_ids[current_id] ";" $9;
    current_id++;
    print $0;
    next;
}
{
    print $0;  # Print other lines unchanged
}' "${TEMP_DIR}/${BASENAME}_prodigal.gff" > "${TEMP_DIR}/${BASENAME}_fixed.gff"

# Step 4: Run PADLOC with the fixed files
echo "[$(date +%T)] Running PADLOC analysis..."
if [ $FORCE -eq 1 ]; then
    FORCE_FLAG="--force"
else
    FORCE_FLAG=""
fi

# Create a temporary output directory for PADLOC
PADLOC_TEMP_DIR="${TEMP_DIR}/padloc_out"
mkdir -p "$PADLOC_TEMP_DIR"

# Run PADLOC with the temporary directory
padloc --faa "${TEMP_DIR}/${BASENAME}_prodigal.faa" --gff "${TEMP_DIR}/${BASENAME}_fixed.gff" --outdir "$PADLOC_TEMP_DIR" --cpu "$CPU" $FORCE_FLAG

# Check if PADLOC was successful
if [ -f "${PADLOC_TEMP_DIR}/${BASENAME}_prodigal_padloc.csv" ]; then
    # Copy the results to the desired output location with the correct name
    cp "${PADLOC_TEMP_DIR}/${BASENAME}_prodigal_padloc.csv" "${OUTPUT_DIR}/${BASENAME}_padloc.csv"
    cp "${PADLOC_TEMP_DIR}/${BASENAME}_prodigal_padloc.gff" "${OUTPUT_DIR}/${BASENAME}_padloc.gff" 2>/dev/null || true
    
    echo "[$(date +%T)] PADLOC analysis completed successfully!"
    echo "Results saved to: ${OUTPUT_DIR}/${BASENAME}_padloc.csv"
    exit 0
else
    echo "[$(date +%T)] Error: PADLOC analysis failed."
    exit 1
fi
