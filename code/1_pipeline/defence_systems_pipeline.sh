#!/bin/bash
# Master pipeline script for Acinetobacter defence system and resistance analysis


echo "===== Acinetobacter Defence System and Resistance Analysis Pipeline ====="
echo "Starting pipeline at $(date)"

# Step 1: Download genomes
echo "Step 1: Downloading genome sequences..."
./code/pipeline/1_download_genomes.sh
if [ $? -ne 0 ]; then
    echo "Error: Genome download failed. Exiting pipeline."
    exit 1
fi

# Step 2: Run DefenseFinder
echo "Step 2: Running DefenseFinder analysis..."
./code/pipeline/2_run_defensefinder.sh
if [ $? -ne 0 ]; then
    echo "Error: DefenseFinder analysis failed. Exiting pipeline."
    exit 1
fi

# Step 3: Run PADLOC
echo "Step 3: Running PADLOC analysis..."
./code/pipeline/3_run_padloc.sh
if [ $? -ne 0 ]; then
    echo "Error: PADLOC analysis failed. Exiting pipeline."
    exit 1
fi

# Step 4: Run ResFinder
echo "Step 4: Running ResFinder analysis..."
./code/pipeline/4_run_resfinder.sh
if [ $? -ne 0 ]; then
    echo "Error: ResFinder analysis failed. Exiting pipeline."
    exit 1
fi

# Step 5: Consolidate results
echo "Step 5: Consolidating all results..."
./code/pipeline/5_consolidate_results.sh
./code/pipeline/6_consolidate_resfinder.sh
if [ $? -ne 0 ]; then
    echo "Error: Result consolidation failed. Exiting pipeline."
    exit 1
fi

# Step 6: Run R analysis
echo "Step 6: Running R analysis..."
Rscript code/analysis/defense_system_analysis.R
if [ $? -ne 0 ]; then
    echo "Error: R analysis failed. Exiting pipeline."
    exit 1
fi

echo "Pipeline completed successfully at $(date)"
