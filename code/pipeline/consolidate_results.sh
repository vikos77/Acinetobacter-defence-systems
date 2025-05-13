#!/bin/bash
# Script to consolidate DefenseFinder and PADLOC results
# Author: Vigneshwaran Muthuraman

# Set directories
DEFENSE_DIR="results/defensefinder"
PADLOC_DIR="results/padloc"
OUTPUT_DIR="results/consolidated"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "Consolidating DefenseFinder results..."

# Initialize a flag to capture header only once
first_df_header_written=false

# Find all defense_finder_systems.tsv files
find "$DEFENSE_DIR" -name "*_defense_finder_systems.tsv" | sort | while read -r FILE; do
    # Extract genome ID from directory path
    DIR_NAME=$(dirname "$FILE")
    GENOME_ID=$(basename "$DIR_NAME")

    if [[ "$first_df_header_written" = false ]]; then
        # Read header from the very first file and write to output
        head -n 1 "$FILE" > "$OUTPUT_DIR/consolidated_defense_systems.tsv"
        first_df_header_written=true
    fi

    # Append the rest of the file (skipping header), prepending genome ID
    tail -n +2 "$FILE" | awk -v genome="$GENOME_ID" '{print genome"\t"$0}' \
        >> "$OUTPUT_DIR/consolidated_defense_systems.tsv"
done

echo "Consolidating PADLOC results..."

# Initialize a flag to capture header only once
first_padloc_header_written=false

# Find all padloc.csv files
find "$PADLOC_DIR" -name "*_padloc.csv" | sort | while read -r FILE; do
    # Extract genome ID from directory path
    DIR_NAME=$(dirname "$FILE")
    GENOME_ID=$(basename "$DIR_NAME")

    if [[ "$first_padloc_header_written" = false ]]; then
        # Read header from the first CSV and write to output
        head -n 1 "$FILE" > "$OUTPUT_DIR/consolidated_padloc_results.csv"
        first_padloc_header_written=true
    fi

    # Append data lines (skipping header), prepending genome ID
    tail -n +2 "$FILE" | awk -v genome="$GENOME_ID" '{OFS="\t";print genome, $0}' \
        >> "$OUTPUT_DIR/consolidated_padloc_results.csv"
done

echo "Results consolidated successfully."
echo "DefenseFinder results: $OUTPUT_DIR/consolidated_defense_systems.tsv"
echo "PADLOC results:     $OUTPUT_DIR/consolidated_padloc_results.csv"

