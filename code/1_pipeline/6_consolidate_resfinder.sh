#!/bin/bash
# Script to consolidate ResFinder results

# Set directories
RESFINDER_DIR="results/resfinder"
OUTPUT_DIR="results/consolidated"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "Consolidating ResFinder results..."

# Gather files into an array 
mapfile -t FILES < <(find "$RESFINDER_DIR" -name "ResFinder_results_tab.txt" | sort)

header_written=false

for FILE in "${FILES[@]}"; do
    GENOME_ID=$(basename "$(dirname "$FILE")")

    # Write header (with Genome_ID) only once
    if ! $header_written; then
        head -n 1 "$FILE" \
          | awk 'BEGIN { OFS="\t" } { print "Genome_ID", $0 }' \
          > "$OUTPUT_DIR/consolidated_resfinder_results.tsv"
        header_written=true
    fi

    # Append data lines, skipping the file’s own header
    tail -n +2 "$FILE" \
      | awk -v genome="$GENOME_ID" 'BEGIN { OFS="\t" } { print genome, $0 }' \
      >> "$OUTPUT_DIR/consolidated_resfinder_results.tsv"
done

echo "ResFinder results consolidated to: $OUTPUT_DIR"
