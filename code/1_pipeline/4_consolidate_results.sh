#!/bin/bash
# Script to consolidate DefenseFinder and PADLOC results

DEFENSE_DIR="results/defensefinder"
PADLOC_DIR="results/padloc"
OUTPUT_DIR="results/consolidated"
mkdir -p "$OUTPUT_DIR"

### Consolidate DefenseFinder ###
echo "Consolidating DefenseFinder results..."

DF_OUT="$OUTPUT_DIR/consolidated_defense_systems.tsv"

# Gather files
mapfile -t DF_FILES < <(find "$DEFENSE_DIR" -name "*_defense_finder_systems.tsv" | sort)

header_written=false
for FILE in "${DF_FILES[@]}"; do
    GENOME_ID=$(basename "$(dirname "$FILE")")

    if ! $header_written; then
        # Prepend Genome_ID to the first header
        head -n 1 "$FILE" \
          | awk 'BEGIN { OFS="\t" } { print "Genome_ID", $0 }' \
          > "$DF_OUT"
        header_written=true
    fi

    # Append data lines, skipping header
    tail -n +2 "$FILE" \
      | awk -v genome="$GENOME_ID" 'BEGIN { OFS="\t" } { print genome, $0 }' \
      >> "$DF_OUT"
done

### Consolidate PADLOC ###
echo "Consolidating PADLOC results..."
PL_OUT="$OUTPUT_DIR/consolidated_padloc_results.tsv"

# Gather files
mapfile -t PL_FILES < <(find "$PADLOC_DIR" -name "*_padloc.csv" | sort)

header_written=false
for FILE in "${PL_FILES[@]}"; do
    GENOME_ID=$(basename "$(dirname "$FILE")")

    if ! $header_written; then
        # Read the comma-delimited header, convert to TSV, prepend Genome_ID
        head -n 1 "$FILE" \
          | sed 's/,/\t/g' \
          | awk 'BEGIN { OFS="\t" } { print "Genome_ID", $0 }' \
          > "$PL_OUT"
        header_written=true
    fi

    # Convert CSV data to TSV, prepend genome ID
    tail -n +2 "$FILE" \
      | sed 's/,/\t/g' \
      | awk -v genome="$GENOME_ID" 'BEGIN { OFS="\t" } { print genome, $0 }' \
      >> "$PL_OUT"
done

echo "Results consolidated successfully."
echo "DefenseFinder results: $DF_OUT"
echo "PADLOC results:       $PL_OUT"

