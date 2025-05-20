# ==============================================================================
# 6_ime_analysis.R
# 
# Analysis of Integrative Mobile Elements (IMEs) in Acinetobacter species and
# their correlation with defense systems
# 
# Author: Vigneshwaran Muthuraman
# ==============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(viridis)

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
ime_blast_file <- "results/ime_analysis/blast_results/IME_proteins_vs_acinetobacter.tblastn"
defensefinder_file <- "results/consolidated/consolidated_defense_systems.tsv"
metadata_file <- "data/metadata/genome_metadata.xlsx"
output_dir <- "results/figures"
ime_output_dir <- "results/ime_analysis/results"

# Create output directories if they don't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
if (!dir.exists(ime_output_dir)) {
  dir.create(ime_output_dir, recursive = TRUE)
}

# ==============================================================================
# Data loading and cleaning
# ==============================================================================

# Load IME BLAST data
ime_blast_data <- read.delim(ime_blast_file, header = FALSE)

# Rename columns
col_names <- c("QueryID", "SubjectID", "PercentIdentity", "AlignmentLength", 
               "Mismatches", "GapOpens", "QueryStart", "QueryEnd", 
               "SubjectStart", "SubjectEnd", "Evalue", "BitScore")
names(ime_blast_data)[1:min(ncol(ime_blast_data), length(col_names))] <- col_names


# Filter hits by quality
filtered_hits <- ime_blast_data %>%
  filter(PercentIdentity >= 80, Evalue <= 1e-6)

# Extract detailed information from query IDs
parsed_hits <- filtered_hits %>%
  mutate(
    # Full query header for reference
    full_header = QueryID,
    
    # Extract ICEberg element ID
    element_id = sub("(ICEberg\\|\\d+).*", "\\1", QueryID),
    
    # Simplify the genome ID 
    genome_id = str_extract(SubjectID, "^\\S+")
  )

# Further parse protein details
parsed_hits <- parsed_hits %>% 
  mutate(
    # Extract protein accession - look for common accession patterns
    protein_accession = case_when(
      # Handle WP_XXXXXX.X format
      grepl("WP_\\d+\\.\\d+", full_header) ~ 
        str_extract(full_header, "WP_\\d+\\.\\d+"),
      # Handle AAX/ATB/etc format (common GenBank/RefSeq)
      grepl("[A-Z]{3}\\d+\\.\\d+", full_header) ~ 
        str_extract(full_header, "[A-Z]{3}\\d+\\.\\d+"),
      # GenBank specific
      grepl("gb\\|([^\\|]+)", full_header) ~ 
        sub(".*gb\\|([^\\|]+)\\|.*", "\\1", full_header),
      # RefSeq specific
      grepl("ref\\|([^\\|]+)", full_header) ~ 
        sub(".*ref\\|([^\\|]+)\\|.*", "\\1", full_header),
      # Default
      TRUE ~ "unknown"
    ),
    
    # Extract protein function
    raw_function = sub(".*\\|(ref|gb|emb)\\|[^\\|]+\\|(.*)\\[.*", "\\2", full_header)
  )

# Clean up protein function and extract more information
parsed_hits <- parsed_hits %>%
  mutate(
    # Clean up protein function - remove pipes and extra spaces
    protein_function = ifelse(raw_function != full_header, 
                              gsub("\\|", " ", raw_function),
                              "unknown"),
    # Extract organism if present
    organism = ifelse(
      grepl("\\[.*\\]", full_header),
      # 1) grab everything inside the brackets
      sub(".*\\[(.*)\\].*", "\\1", full_header) %>%
        # 2) replace internal '|' with spaces
        gsub("\\|", " ", .) %>%
        # 3) trim any leading/trailing whitespace
        trimws(),
      "unknown"
    )
  )

# Save parsed hits for reference
saveRDS(parsed_hits, file.path(ime_output_dir, "parsed_hits.rds"))

# Group by element_id and protein_function to count prevalence
protein_prevalence <- parsed_hits %>%
  group_by(element_id, protein_function, protein_accession) %>%
  summarize(
    genome_count = n_distinct(genome_id),
    total_hits = n(),
    mean_identity = mean(PercentIdentity),
    best_evalue = min(Evalue),
    organism = first(organism),
    .groups = "drop"
  ) %>%
  arrange(desc(genome_count), desc(total_hits))

# Count element prevalence for context
element_prevalence <- parsed_hits %>%
  group_by(element_id) %>%
  summarize(
    genome_count = n_distinct(genome_id),
    protein_count = n_distinct(protein_function),
    .groups = "drop"
  ) %>%
  arrange(desc(genome_count))

# Save summaries
write_csv(protein_prevalence, file.path(ime_output_dir, "protein_prevalence.csv"))
write_csv(element_prevalence, file.path(ime_output_dir, "element_prevalence.csv"))

# Load DefenseFinder data
defense_df <- read_tsv(defensefinder_file, show_col_types = FALSE)

# Clean defense data - ensure unique defense systems per genome
defense_clean <- defense_df %>%
  distinct(Genome_ID, type)

# ==============================================================================
# Panel A: Distribution of IME Counts per Genome
# ==============================================================================

# Count IMEs per genome
ime_counts <- parsed_hits %>%
  distinct(genome_id, element_id) %>%
  count(genome_id, name = "ime_count")

# Create histogram
panel_a <- ggplot(ime_counts, aes(x = ime_count)) +
  geom_histogram(binwidth = 1, fill = "#FF7F00", color = "black", alpha = 0.8) +
  geom_vline(aes(xintercept = mean(ime_count, na.rm = TRUE)), 
             linetype = "dashed", color = "darkorange3", size = 1) +
  annotate("text", x = mean(ime_counts$ime_count) + 1, 
           y = max(table(ime_counts$ime_count)) * 0.9, 
           label = paste("Mean:", round(mean(ime_counts$ime_count, na.rm = TRUE), 1)),
           color = "darkorange3", fontface = "bold") +
  labs(
    title = "A. Distribution of IME Element Counts",
    x = "Number of IME Elements per Genome",
    y = "Number of Genomes"
  ) +
  custom_theme

# ==============================================================================
# Panel B: Most Prevalent IME Elements
# ==============================================================================

# Create bar chart of top IME elements
panel_b <- element_prevalence %>%
  head(20) %>%
  mutate(
    # Format element_id for better readability
    display_id = str_extract(element_id, "\\d+$"),
    display_id = paste0("ICEberg-", display_id),
    display_id = fct_reorder(display_id, genome_count)
  ) %>%
  ggplot(aes(x = display_id, y = genome_count)) +
  geom_bar(stat = "identity", fill = "#FF7F00", alpha = 0.8) +
  geom_text(aes(label = genome_count), hjust = -0.2, size = 3.5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "B. Top 20 IME Elements",
    x = NULL,
    y = "Number of Genomes"
  ) +
  custom_theme +
  theme(axis.text.y = element_text(face = "italic"))

# ==============================================================================
# Panel C: Defense vs IME Count Correlation
# ==============================================================================

# Count defense systems per genome
defense_counts <- defense_clean %>%
  count(Genome_ID, name = "defense_count")

# Combine defense and IME counts
combined_counts <- ime_counts %>%
  rename(Genome_ID = genome_id) %>%
  full_join(defense_counts, by = "Genome_ID") %>%
  # Fill NAs with zeros
  mutate(
    defense_count = replace_na(defense_count, 0),
    ime_count = replace_na(ime_count, 0)
  )

# Calculate correlation
corr_test <- cor.test(combined_counts$defense_count, 
                      combined_counts$ime_count, 
                      method = "spearman")

# Create scatter plot
panel_c <- ggplot(combined_counts, aes(x = defense_count, y = ime_count)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "darkorange3") +
  labs(
    title = "C. Correlation Between Defense Systems and IMEs",
    subtitle = sprintf("Spearman's ρ = %.2f, p = %.4f", 
                       corr_test$estimate, 
                       corr_test$p.value),
    x = "Number of Defense Systems",
    y = "Number of IME Elements"
  ) +
  custom_theme

# ==============================================================================
# Panel D: Correlation Between Defense Systems and IME Proteins
# ==============================================================================

# Use the correlation analysis function from your script
analyze_defense_protein_correlations <- function(defense_df, ime_protein_data, output_dir, n_top = 10) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Step 1: Get the top n defense systems by prevalence
  defense_counts <- defense_df %>%
    group_by(type) %>% 
    count(type, sort = TRUE) %>%
    head(n_top)
  
  cat("Top", n_top, "defense systems by prevalence:\n")
  print(defense_counts)
  
  top_defense <- defense_counts$type
  
  # Step 2: Get the top n protein functions by genome count
  top_proteins <- ime_protein_data %>%
    arrange(desc(genome_count)) %>%
    distinct(element_id, protein_function, .keep_all = TRUE) %>%
    head(10)
  
  cat("\nTop 10 IME protein functions by genome prevalence:\n")
  print(top_proteins[, c("element_id", "protein_function", "genome_count")])
  
  # Save the top protein info to CSV
  write.csv(top_proteins, file.path(output_dir, "top_protein_functions.csv"), row.names = FALSE)
  
  # Step 3: Create binary matrices for defense systems and protein functions
  
  # For defense systems - get unique genome-defense pairs
  defense_binary <- defense_df %>%
    filter(type %in% top_defense) %>%
    mutate(present = 1) %>%
    distinct(Genome_ID, type, present) %>%
    pivot_wider(
      names_from = type,
      values_from = present,
      values_fill = 0
    )
  
  # For protein functions - create a binary matrix from parsed_hits data
  protein_binary_data <- data.frame()
  
  # For each top protein, find all genomes with that protein
  for (i in 1:nrow(top_proteins)) {
    current_element <- top_proteins$element_id[i]
    current_function <- top_proteins$protein_function[i]
    current_accession <- top_proteins$protein_accession[i]
    
    # Find all genomes with this protein
    matching_genomes <- parsed_hits %>%
      filter(
        element_id == current_element & 
          grepl(current_function, protein_function, fixed = TRUE) &
          (protein_accession == current_accession | is.na(current_accession))
      ) %>%
      pull(genome_id) %>%
      unique()
    
    # Clean up the function name for use in column names
    clean_function <- gsub("[^a-zA-Z0-9]", "_", substr(current_function, 1, 30))
    
    # Create column name with element ID and function
    col_name <- paste0(current_element, "_", clean_function)
    
    # Add to our dataframe
    if (length(matching_genomes) > 0) {
      new_rows <- data.frame(
        genome_id = matching_genomes,
        col_name = rep(1, length(matching_genomes))
      )
      names(new_rows)[2] <- col_name
      
      if (nrow(protein_binary_data) == 0) {
        protein_binary_data <- new_rows
      } else {
        protein_binary_data <- full_join(protein_binary_data, new_rows, by = "genome_id")
      }
    }
  }
  
  # Fill NA values with 0
  protein_binary_data[is.na(protein_binary_data)] <- 0
  
  # Step 4: Merge the two datasets
  # First make sure we're using consistent genome ID column names
  names(defense_binary)[1] <- "genome_id"
  
  # Merge the datasets
  merged_data <- defense_binary %>%
    full_join(protein_binary_data, by = "genome_id") %>%
    # Replace NA values with 0
    mutate(across(everything(), ~ifelse(is.na(.), 0, .)))
  
  # Remove the ID column and ensure all data is numeric
  correlation_data <- merged_data %>%
    select(-genome_id) %>%
    mutate(across(everything(), as.numeric))
  
  # Check for any columns with zero variance
  zero_var <- apply(correlation_data, 2, function(x) length(unique(x)) <= 1)
  if (any(zero_var)) {
    cat("Removing", sum(zero_var), "columns with zero variance\n")
    correlation_data <- correlation_data[, !zero_var]
  }
  
  # Define defense and protein columns
  defense_cols <- intersect(names(correlation_data), top_defense)
  protein_cols <- setdiff(names(correlation_data), defense_cols)
  
  # Step 5: Perform Fisher's exact tests
  cat("Performing Fisher's exact tests...\n")
  
  # Initialize results dataframe
  fisher_results <- data.frame(
    Defense_System = character(),
    Protein_Function = character(),
    Odds_Ratio = numeric(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # For each defense system-protein pair
  for (defense in defense_cols) {
    for (protein in protein_cols) {
      # Create 2x2 contingency table
      n11 <- sum(correlation_data[, defense] == 1 & correlation_data[, protein] == 1)  # Both present
      n10 <- sum(correlation_data[, defense] == 1 & correlation_data[, protein] == 0)  # Defense present, protein absent
      n01 <- sum(correlation_data[, defense] == 0 & correlation_data[, protein] == 1)  # Defense absent, protein present
      n00 <- sum(correlation_data[, defense] == 0 & correlation_data[, protein] == 0)  # Both absent
      
      # Create contingency table
      contingency_table <- matrix(c(n11, n01, n10, n00), nrow = 2, byrow = TRUE)
      
      # Perform Fisher's exact test
      tryCatch({
        test_result <- fisher.test(contingency_table)
        
        # Create a more readable function name for display
        readable_function <- sub(".*_", "", protein)
        
        # Add results to dataframe
        fisher_results <- rbind(fisher_results, data.frame(
          Defense_System = defense,
          Protein_Function = protein,
          Readable_Function = readable_function,
          Odds_Ratio = test_result$estimate,
          P_Value = test_result$p.value,
          stringsAsFactors = FALSE
        ))
      }, error = function(e) {
        cat("Error testing", defense, "vs", protein, ":", conditionMessage(e), "\n")
      })
    }
  }
  
  # Apply FDR correction
  if (nrow(fisher_results) > 0) {
    fisher_results$P_Adjusted <- p.adjust(fisher_results$P_Value, method = "BH")
    
    # Add log2 odds ratio for better visualization
    fisher_results$log_odds <- log2(fisher_results$Odds_Ratio)
    
    # Handle infinite values
    fisher_results$log_odds[is.infinite(fisher_results$log_odds) & fisher_results$log_odds > 0] <- 16
    fisher_results$log_odds[is.infinite(fisher_results$log_odds) & fisher_results$log_odds < 0] <- -16
    
    # Add significance levels
    fisher_results$Significance <- ""
    fisher_results$Significance[fisher_results$P_Adjusted < 0.05] <- "*"
    fisher_results$Significance[fisher_results$P_Adjusted < 0.01] <- "**"
    fisher_results$Significance[fisher_results$P_Adjusted < 0.001] <- "***"
    
    # Save Fisher's test results
    write.csv(fisher_results, file.path(output_dir, "defense_protein_fisher_tests.csv"), row.names = FALSE)
  } else {
    cat("No valid Fisher's test results obtained\n")
    return(NULL)
  }
  
  # Step 6: Create heatmap visualization for odds ratios
  if (requireNamespace("pheatmap", quietly = TRUE) && nrow(fisher_results) > 0) {
    # Create odds ratio matrix for visualization
    odds_matrix <- matrix(NA, nrow = length(defense_cols), ncol = length(protein_cols))
    rownames(odds_matrix) <- defense_cols
    colnames(odds_matrix) <- protein_cols
    
    # Fill the matrix with log2 odds ratios
    for (i in 1:nrow(fisher_results)) {
      defense <- fisher_results$Defense_System[i]
      protein <- fisher_results$Protein_Function[i]
      odds_matrix[defense, protein] <- fisher_results$log_odds[i]
    }
    
    # Create significance indicators matrix
    labels_matrix <- matrix("", nrow = length(defense_cols), ncol = length(protein_cols))
    rownames(labels_matrix) <- defense_cols
    colnames(labels_matrix) <- protein_cols
    
    for (i in 1:nrow(fisher_results)) {
      defense <- fisher_results$Defense_System[i]
      protein <- fisher_results$Protein_Function[i]
      labels_matrix[defense, protein] <- fisher_results$Significance[i]
    }
    
    # Create readable column names by extracting element ID and protein function
    readable_cols <- sapply(protein_cols, function(p) {
      # Extract ICEberg ID without the "ICEberg|" prefix
      iceberg_id <- sub("^ICEberg\\|(\\d+)_.*", "\\1", p)
      
      # Extract protein function (everything after the element_id and underscore)
      func_part <- sub("^ICEberg\\|\\d+_(.*)", "\\1", p)
      
      # Clean up function part (replace underscores with spaces)
      func_part <- gsub("_", " ", func_part)
      
      # Truncate if too long
      if (nchar(func_part) > 18) {
        func_part <- paste0(substr(func_part, 1, 15), "...")
      }
      
      # Combine ID and function
      return(paste0(iceberg_id, "-", func_part))
    })
    
    # Apply the readable column names
    colnames(odds_matrix) <- readable_cols
    colnames(labels_matrix) <- readable_cols
    
    # Ensure we close any open graphics devices
    while (!is.null(dev.list())) {
      dev.off()
    }
    
    # Create pheatmap with better formatting
    pheatmap::pheatmap(
      odds_matrix,
      display_numbers = labels_matrix,
      color = colorRampPalette(c("#E41A1C", "white", "blue"))(100),
      breaks = seq(-8, 8, length.out = 101),
      main = "D. Association Between Defense Systems and IME Proteins",
      filename = file.path(output_dir, "defense_ime_heatmap.png"),
      width = 12,  
      height = 8,
      fontsize = 10,
      fontsize_number = 12,
      fontsize_col = 8,  
      number_color = "black",
      border_color = "white",
      cellwidth = 35,
      cellheight = 35,
      na_col = "grey90",
      angle_col = 45  # Angle the column names for better readability
    )
    
    # Save PDF version
    pheatmap::pheatmap(
      odds_matrix,
      display_numbers = labels_matrix,
      color = colorRampPalette(c("#E41A1C", "white", "blue"))(100),
      breaks = seq(-8, 8, length.out = 101),
      main = "D. Association Between Defense Systems and IME Proteins",
      filename = file.path(output_dir, "defense_ime_heatmap.pdf"),
      width = 12,  
      height = 8,
      fontsize = 10,
      fontsize_number = 12,
      fontsize_col = 8,  
      number_color = "black",
      border_color = "white",
      cellwidth = 35,
      cellheight = 35,
      na_col = "grey90",
      angle_col = 45
    )
    
    # Save the column name mapping for reference
    name_mapping <- data.frame(
      original_name = protein_cols,
      readable_name = readable_cols
    )
    write.csv(name_mapping, file.path(output_dir, "protein_name_mapping.csv"), row.names = FALSE)
  }
  
  # Return results
  return(list(
    fisher_results = fisher_results,
    top_defense = defense_cols,
    top_proteins = protein_cols,
    correlation_data = correlation_data
  ))
}

# Run correlation analysis
correlation_results <- analyze_defense_protein_correlations(
  defense_df = defense_clean,
  ime_protein_data = protein_prevalence,
  output_dir = ime_output_dir,
  n_top = 10
)

# Create panel D from the visualization output
panel_d <- grid::rasterGrob(
  png::readPNG(file.path(ime_output_dir, "defense_ime_heatmap.png")),
  interpolate = TRUE
)

# ==============================================================================
# Combine panels and save figure
# ==============================================================================

# Save individual panels
ggsave(file.path(output_dir, "ime_panel_a.png"), panel_a, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "ime_panel_b.png"), panel_b, width = 6, height = 6, dpi = 300)
ggsave(file.path(output_dir, "ime_panel_c.png"), panel_c, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "ime_panel_d.png"), panel_d, width = 6, height = 5, dpi = 300)
# Panel D already saved by the pheatmap function

# Combine all panels into a single figure
combined_figure <- grid.arrange(
  panel_a, panel_b,
  panel_c, panel_d,
  ncol = 2,
  nrow = 2
)

# Save combined figure
ggsave(file.path(output_dir, "Figure11_ime_analysis.png"), combined_figure, width = 12, height = 10, dpi = 300)
ggsave(file.path(output_dir, "Figure11_ime_analysis.pdf"), combined_figure, width = 12, height = 10)

