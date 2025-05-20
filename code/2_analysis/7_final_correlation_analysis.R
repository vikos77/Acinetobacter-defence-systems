# ==============================================================================
# 7_final_correlation_analysis.R
# 
# Final correlation analysis between all genomic elements:
# Defense systems, Antibiotic Resistance Genes (ARGs), Anti-defense systems, 
# and Integrative Mobile Elements (IMEs)
# 
# Author: Vigneshwaran Muthuraman
# ==============================================================================

# Load required libraries
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(png)

# Set theme for consistent visualization
theme_set(theme_bw(base_size = 12))
custom_theme <- theme(
  axis.title = element_text(face = "bold"),
  axis.text = element_text(size = 10),
  plot.title = element_text(size = 12, face = "bold"),
  legend.title = element_text(size = 10, face = "bold"),
  legend.text = element_text(size = 9),
  panel.grid.minor = element_blank()
)

# Set file paths
defensefinder_file <- "results/consolidated/consolidated_defense_systems.tsv"
resfinder_file <- "results/consolidated/consolidated_resfinder_results.tsv"
antidefense_file <- "results/consolidated/consolidated_antidefense_systems.tsv"
ime_blast_file <- "results/ime_analysis/blast_results/IME_proteins_vs_acinetobacter.tblastn"
output_dir <- "results/figures"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# Enhanced function to create correlation heatmap with proper data handling
# ==============================================================================

create_correlation_heatmap <- function(defense_df, arg_df, anti_defense_df, ime_blast_file, 
                                       output_file = "correlation_heatmap.png") {
  

  # 1. DEFENSE SYSTEMS - count distinct defense types per genome

  defense_counts <- defense_df %>%
    select(Genome_ID, type) %>%
    # Remove duplicates: same genome with same defense system type
    distinct() %>%
    group_by(Genome_ID) %>%
    summarize(defense_count = n(), .groups = "drop")
  
  cat("Defense systems: Found", nrow(defense_counts), "genomes with mean count =", 
      round(mean(defense_counts$defense_count), 2), "\n")
  
  # 2. ARGS - count distinct resistance genes per genome

  if ("arg_count" %in% colnames(arg_df)) {
    # If already processed
    arg_counts <- arg_df %>% select(Genome_ID, arg_count)
  } else if ("Resistance gene" %in% colnames(arg_df)) {
    arg_counts <- arg_df %>%
      select(Genome_ID, `Resistance gene`) %>%
      # Remove duplicates: same genome with same resistance gene
      distinct() %>%
      group_by(Genome_ID) %>%
      summarize(arg_count = n(), .groups = "drop")
  } else {
    stop("Cannot identify resistance gene column in arg_df")
  }
  
  cat("ARGs: Found", nrow(arg_counts), "genomes with mean count =", 
      round(mean(arg_counts$arg_count), 2), "\n")
  
  # 3. ANTI-DEFENSE SYSTEMS - count distinct anti-defense types per genome
  message("Processing anti-defense systems...")
  if ("type" %in% colnames(anti_defense_df)) {
    anti_defense_counts <- anti_defense_df %>%
      select(Genome_ID, type) %>%
      # Remove duplicates: same genome with same anti-defense system type
      distinct() %>%
      group_by(Genome_ID) %>%
      summarize(anti_defense_count = n(), .groups = "drop")
  } else {
    stop("Cannot find 'type' column in anti_defense_df")
  }
  
  cat("Anti-defense systems: Found", nrow(anti_defense_counts), "genomes with mean count =", 
      round(mean(anti_defense_counts$anti_defense_count), 2), "\n")
  
  # 4. IME - count distinct IME elements per genome from tblastn results
  message("Processing integrative mobile elements...")
  # Read tblastn results (protein query vs nucleotide database)
  blast_cols <- c("QueryID", "SubjectID", "PercentIdentity", "AlignmentLength", 
                  "Mismatches", "GapOpens", "QueryStart", "QueryEnd", 
                  "SubjectStart", "SubjectEnd", "Evalue", "BitScore")
  
  # Check if file exists
  if (!file.exists(ime_blast_file)) {
    warning("IME BLAST file not found. Setting IME counts to zero for all genomes.")
    # Create empty IME data
    all_genomes <- unique(c(defense_counts$Genome_ID, arg_counts$Genome_ID, anti_defense_counts$Genome_ID))
    ime_data <- data.frame(
      Genome_ID = all_genomes,
      ime_count = 0
    )
  } else {
    # Read and process BLAST data
    blast_data <- tryCatch({
      read.delim(ime_blast_file, header = FALSE, stringsAsFactors = FALSE)
    }, error = function(e) {
      warning("Error reading IME BLAST file: ", e$message)
      return(data.frame())
    })
    
    if (nrow(blast_data) > 0) {
      # Assign column names based on available columns
      n_cols <- min(ncol(blast_data), length(blast_cols))
      colnames(blast_data)[1:n_cols] <- blast_cols[1:n_cols]
      
      # Process tblastn results similar to your IME analysis
      ime_data <- blast_data %>%
        # Filter for high quality hits
        filter(PercentIdentity >= 80, Evalue <= 1e-6) %>%
        mutate(
          # Extract genome ID from SubjectID
          Genome_ID = str_extract(SubjectID, "^\\S+"),
          
          # Extract ICEberg element ID from QueryID
          # This extracts the ICEberg ID pattern: ICEberg|number
          element_id = sub("(ICEberg\\|\\d+).*", "\\1", QueryID)
        ) %>%
        # Remove hits where we couldn't extract an element_id
        filter(!is.na(element_id) & element_id != "" & element_id != QueryID) %>%
        # Remove duplicates: same genome with same IME element
        select(Genome_ID, element_id) %>%
        distinct() %>%
        # Count distinct IME elements per genome
        group_by(Genome_ID) %>%
        summarize(ime_count = n(), .groups = "drop")
    } else {
      # Create empty IME data if BLAST file is empty
      all_genomes <- unique(c(defense_counts$Genome_ID, arg_counts$Genome_ID, anti_defense_counts$Genome_ID))
      ime_data <- data.frame(
        Genome_ID = all_genomes,
        ime_count = 0
      )
    }
  }
  
  cat("IMEs: Found", nrow(ime_data), "genomes with mean count =", 
      round(mean(ime_data$ime_count), 2), "\n")
  
  # 5. MERGE ALL COUNTS - ensure all are numeric
  combined_data <- defense_counts %>%
    full_join(arg_counts, by = "Genome_ID") %>%
    full_join(anti_defense_counts, by = "Genome_ID") %>%
    full_join(ime_data, by = "Genome_ID") %>%
    # Replace NA with 0 and ensure all counts are numeric
    mutate(
      defense_count = as.numeric(ifelse(is.na(defense_count), 0, defense_count)),
      arg_count = as.numeric(ifelse(is.na(arg_count), 0, arg_count)),
      anti_defense_count = as.numeric(ifelse(is.na(anti_defense_count), 0, anti_defense_count)),
      ime_count = as.numeric(ifelse(is.na(ime_count), 0, ime_count))
    )
  
  cat("\nCombined data: Total", nrow(combined_data), "genomes\n")
  
  # Show summary statistics
  summary_stats <- combined_data %>%
    summarize(across(ends_with("_count"), 
                     list(mean = ~mean(.), min = ~min(.), max = ~max(.), sd = ~sd(.)), 
                     .names = "{.col}_{.fn}"))
  
  cat("\nSummary Statistics:\n")
  print(summary_stats)
  
  # 6. CALCULATE CORRELATION MATRIX
  # Convert to matrix and ensure numeric
  cor_data <- as.matrix(combined_data[, c("defense_count", "arg_count", "anti_defense_count", "ime_count")])
  colnames(cor_data) <- c("Defence", "ARG", "Anti-Defence", "IME")
  
  # Verify data is numeric
  cat("\nVerifying data types:\n")
  print(sapply(data.frame(cor_data), class))
  
  # Check for zero variance columns
  zero_var_cols <- which(apply(cor_data, 2, var) == 0)
  if (length(zero_var_cols) > 0) {
    cat("Warning: Columns with zero variance:", colnames(cor_data)[zero_var_cols], "\n")
    # Don't remove them, but note them
  }
  
  # Calculate correlation matrix using Spearman correlation
  cor_matrix <- cor(cor_data, method = "spearman", use = "pairwise.complete.obs")
  
  # Print correlation matrix
  cat("\nCorrelation Matrix:\n")
  print(round(cor_matrix, 3))
  
  # 7. CALCULATE P-VALUES FOR CORRELATIONS
  n <- nrow(cor_data)
  p_values <- matrix(NA, nrow = ncol(cor_data), ncol = ncol(cor_data))
  rownames(p_values) <- colnames(p_values) <- colnames(cor_data)
  
  for (i in 1:ncol(cor_data)) {
    for (j in 1:ncol(cor_data)) {
      if (i != j) {
        # Extract vectors properly and ensure they're numeric
        x <- as.numeric(cor_data[,i])
        y <- as.numeric(cor_data[,j])
        test_result <- cor.test(x, y, method = "spearman")
        p_values[i,j] <- test_result$p.value
      } else {
        p_values[i,j] <- 1
      }
    }
  }
  
  cat("\nP-values:\n")
  print(round(p_values, 4))
  
  # 8. CREATE SIGNIFICANCE ANNOTATIONS
  sig_matrix <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  rownames(sig_matrix) <- rownames(cor_matrix)
  colnames(sig_matrix) <- colnames(cor_matrix)
  
  for (i in 1:nrow(p_values)) {
    for (j in 1:ncol(p_values)) {
      if (i != j) {
        if (p_values[i,j] < 0.001) {
          sig_matrix[i,j] <- "***"
        } else if (p_values[i,j] < 0.01) {
          sig_matrix[i,j] <- "**"
        } else if (p_values[i,j] < 0.05) {
          sig_matrix[i,j] <- "*"
        }
      }
    }
  }
  
  # 9. CREATE THE HEATMAP
  # Ensure all graphics devices are closed
  while (!is.null(dev.list())) {
    dev.off()
  }
  
  # Create custom labels that include correlation values and significance
  custom_labels <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  for (i in 1:nrow(cor_matrix)) {
    for (j in 1:ncol(cor_matrix)) {
      if (i == j) {
        custom_labels[i,j] <- sprintf("%.2f", cor_matrix[i,j])
      } else {
        custom_labels[i,j] <- paste0(sprintf("%.2f", cor_matrix[i,j]), "\n", sig_matrix[i,j])
      }
    }
  }
  
  # Create PNG version
  png(output_file, width = 10*300, height = 10*300, res = 300)
  
  pheatmap(
    cor_matrix,
    display_numbers = custom_labels,
    number_color = "black",
    fontsize_number = 14,
    color = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", 
                               "#FFFFFF", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    angle_col = 45,
    main = "D. Correlation Matrix Between Genomic Elements",
    fontsize = 16,
    fontsize_row = 14,
    fontsize_col = 14,
    cellwidth = 80,
    cellheight = 80,
    border_color = "white"
  )
  
  dev.off()
  
  # Create PDF version
  pdf(sub("\\.png$", ".pdf", output_file), width = 10, height = 10)
  
  pheatmap(
    cor_matrix,
    display_numbers = custom_labels,
    number_color = "black",
    fontsize_number = 14,
    color = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", 
                               "#FFFFFF", "#F4A582", "#D6604D", "#B2182B", "#67001F"))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    angle_col = 45,
    main = "Correlation Matrix Between Genomic Elements",
    fontsize = 16,
    fontsize_row = 14,
    fontsize_col = 14,
    cellwidth = 80,
    cellheight = 80,
    border_color = "white"
  )
  
  dev.off()
  
  cat("\nHeatmap saved to:", output_file, "\n")
  
  # 10. SAVE DETAILED RESULTS
  write_csv(as.data.frame(cor_matrix), sub("\\.png$", "_correlations.csv", output_file))
  write_csv(as.data.frame(p_values), sub("\\.png$", "_pvalues.csv", output_file))
  write_csv(combined_data, sub("\\.png$", "_counts_data.csv", output_file))
  
  return(list(correlations = cor_matrix, p_values = p_values, data = combined_data))
}

# ==============================================================================
# Load all data files
# ==============================================================================
# Load DefenseFinder data
defense_df <- read_tsv(defensefinder_file, show_col_types = FALSE)

# Load ResFinder data
resfinder_df <- read_tsv(resfinder_file, show_col_types = FALSE)

# Load AntiDefenseFinder data
antidefense_df <- read_tsv(antidefense_file, show_col_types = FALSE)

# ==============================================================================
# Run the correlation analysis
# ==============================================================================
results <- create_correlation_heatmap(
  defense_df = defense_df,
  arg_df = resfinder_df,
  anti_defense_df = antidefense_df,
  ime_blast_file = ime_blast_file,
  output_file = file.path(output_dir, "final_correlation_matrix.png")
)

# ==============================================================================
# Display and summarize results
# ==============================================================================


cat("\nFinal Correlation Matrix:\n")
print(round(results$correlations, 3))

cat("\nSignificant Correlations (p < 0.05):\n")
significant_mask <- results$p_values < 0.05
diag(significant_mask) <- FALSE  # Remove diagonal

sig_pairs <- which(significant_mask, arr.ind = TRUE)
if (nrow(sig_pairs) > 0) {
  for (i in 1:nrow(sig_pairs)) {
    row_idx <- sig_pairs[i,1]
    col_idx <- sig_pairs[i,2]
    
    # Determine significance level
    p_val <- results$p_values[row_idx, col_idx]
    if (p_val < 0.001) {
      sig_level <- "***"
    } else if (p_val < 0.01) {
      sig_level <- "**"
    } else if (p_val < 0.05) {
      sig_level <- "*"
    } else {
      sig_level <- ""
    }
    
    cat(sprintf("%s vs %s: ρ = %.3f, p = %.4f %s\n", 
                rownames(results$correlations)[row_idx],
                colnames(results$correlations)[col_idx],
                results$correlations[row_idx, col_idx],
                results$p_values[row_idx, col_idx],
                sig_level))
  }
} else {
  cat("No significant correlations found\n")
}