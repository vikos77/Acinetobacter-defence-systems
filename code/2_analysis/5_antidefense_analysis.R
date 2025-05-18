# ==============================================================================
# 5_antidefense_analysis.R
# 
# Analysis of anti-defense systems in Acinetobacter species and
# their correlation with defense systems
# 
# Author: Vigneshwaran Muthuraman
# ==============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(reshape2)

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
antidefense_file <- "results/consolidated/consolidated_antidefense_systems.tsv"
defensefinder_file <- "results/consolidated/consolidated_defense_systems.tsv"
output_dir <- "results/figures"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# Load and prepare data
# ==============================================================================

# Load AntiDefenseFinder data
message("Loading AntiDefenseFinder results...")
anti_defense_df <- read_tsv(antidefense_file, show_col_types = FALSE)

# Load DefenseFinder data for correlation analysis
message("Loading DefenseFinder results...")
defense_df <- read_tsv(defensefinder_file, show_col_types = FALSE)

# Clean DefenseFinder data - ensure unique defense systems per genome
defense_clean <- defense_df %>%
  distinct(Genome_ID, type)

# Count total genomes in dataset
total_genomes <- n_distinct(c(anti_defense_clean$Genome_ID, defense_clean$Genome_ID))
message(paste("Total unique genomes in dataset:", total_genomes))

# ==============================================================================
# Panel A: Distribution of Anti-Defense System Counts per Genome
# ==============================================================================

# Count anti-defense systems per genome
antidefense_counts <- anti_defense_clean %>%
  count(Genome_ID, name = "antidefense_count") %>%
  # Add genomes with zero anti-defense systems, if any
  right_join(data.frame(Genome_ID = unique(c(anti_defense_clean$Genome_ID, defense_clean$Genome_ID))), 
             by = "Genome_ID") %>%
  mutate(antidefense_count = ifelse(is.na(antidefense_count), 0, antidefense_count))

# Create histogram of anti-defense counts
panel_a <- ggplot(antidefense_counts, aes(x = antidefense_count)) +
  geom_histogram(binwidth = 1, fill = "#4DAF4A", color = "black", alpha = 0.8) +
  geom_vline(aes(xintercept = mean(antidefense_count, na.rm = TRUE)), 
             linetype = "dashed", color = "darkgreen", size = 1) +
  annotate("text", x = mean(antidefense_counts$antidefense_count) + 1, 
           y = max(table(antidefense_counts$antidefense_count)) * 0.9, 
           label = paste("Mean:", round(mean(antidefense_counts$antidefense_count, na.rm = TRUE), 1)),
           color = "darkgreen", fontface = "bold") +
  labs(
    title = "A. Distribution of Anti-Defense System Counts",
    x = "Number of Anti-Defense Systems per Genome",
    y = "Number of Genomes"
  ) +
  custom_theme

# ==============================================================================
# Panel B: Prevalence of Different Anti-Defense System Types
# ==============================================================================

# Count and calculate prevalence of each anti-defense system
antidefense_prevalence <- anti_defense_clean %>%
  count(type, name = "count") %>%
  arrange(desc(count)) %>%
  mutate(
    percentage = count / total_genomes * 100,
    type = fct_reorder(type, count, .desc = TRUE)
  )

# Create bar chart
panel_b <- ggplot(antidefense_prevalence, aes(x = type, y = count)) +
  geom_bar(stat = "identity", fill = "#4DAF4A", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "B. Anti-Defense Systems Prevalence",
    x = NULL,
    y = "Number of Genomes"
  ) +
  custom_theme +
  theme(axis.text.y = element_text(face = "italic", size = 9))

# ==============================================================================
# Panel C: Correlation with Defense System Counts
# ==============================================================================

# Count defense systems per genome
defense_counts <- defense_clean %>%
  count(Genome_ID, name = "defense_count")

# Join with anti-defense counts
combined_counts <- antidefense_counts %>%
  left_join(defense_counts, by = "Genome_ID") %>%
  # Fill missing values with zero
  mutate(defense_count = ifelse(is.na(defense_count), 0, defense_count))

# Calculate correlation
correlation_test <- cor.test(combined_counts$defense_count, 
                             combined_counts$antidefense_count, 
                             method = "spearman")

# Create scatter plot with smoothed trend line
panel_c <- ggplot(combined_counts, aes(x = defense_count, y = antidefense_count)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "darkred") +
  labs(
    title = "C. Defense vs Anti-Defense Systems",
    subtitle = sprintf("Spearman's ρ = %.2f, p = %.4f", 
                       correlation_test$estimate, 
                       correlation_test$p.value),
    x = "Number of Defense Systems",
    y = "Number of Anti-Defense Systems"
  ) +
  custom_theme

# ==============================================================================
# Panel D: Defense-Antidefense Associations Heatmap
# ==============================================================================

# Get top defense systems for analysis
top_defense_systems <- defense_clean %>%
  count(type, sort = TRUE) %>%
  head(10) %>%
  pull(type)

# Get all anti-defense systems for analysis
antidefense_systems <- unique(anti_defense_clean$type)

# Function to perform Fisher's exact test for all pairs
perform_fisher_tests <- function(defense_df, anti_defense_df) {
  # Create a dataframe to store results
  fisher_results <- data.frame(
    defense_system = character(),
    antidefense_system = character(),
    odds_ratio = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get all unique genomes
  all_genomes <- unique(c(defense_df$Genome_ID, anti_defense_df$Genome_ID))
  total_genomes <- length(all_genomes)
  
  # Loop through each defense-antidefense pair
  for (def_sys in top_defense_systems) {
    # Get genomes with this defense system
    def_genomes <- unique(defense_df$Genome_ID[defense_df$type == def_sys])
    
    for (antidef_sys in antidefense_systems) {
      # Get genomes with this antidefense system
      antidef_genomes <- unique(anti_defense_df$Genome_ID[anti_defense_df$type == antidef_sys])
      
      # Create 2x2 contingency table
      n11 <- length(intersect(def_genomes, antidef_genomes))  # Both present
      n10 <- length(def_genomes) - n11  # Defense present, antidefense absent
      n01 <- length(antidef_genomes) - n11  # Defense absent, antidefense present
      n00 <- total_genomes - n11 - n10 - n01  # Both absent
      
      # Create the contingency table
      contingency_table <- matrix(c(n11, n01, n10, n00), nrow = 2, byrow = TRUE)
      
      # Perform Fisher's exact test
      tryCatch({
        test_result <- fisher.test(contingency_table)
        
        # Add results to the dataframe
        fisher_results <- rbind(fisher_results, data.frame(
          defense_system = def_sys,
          antidefense_system = antidef_sys,
          odds_ratio = test_result$estimate,
          p_value = test_result$p.value,
          stringsAsFactors = FALSE
        ))
      }, error = function(e) {
        message("Error in Fisher's test for ", def_sys, " vs ", antidef_sys, ": ", conditionMessage(e))
      })
    }
  }
  
  # Apply FDR correction
  fisher_results$p_adjusted <- p.adjust(fisher_results$p_value, method = "BH")
  
  # Add log2 odds ratio for better visualization
  fisher_results$log_odds <- log2(fisher_results$odds_ratio)
  
  # Handle infinite values
  fisher_results$log_odds[is.infinite(fisher_results$log_odds) & fisher_results$log_odds > 0] <- 16
  fisher_results$log_odds[is.infinite(fisher_results$log_odds) & fisher_results$log_odds < 0] <- -16
  
  # Add association direction
  fisher_results$association <- ifelse(fisher_results$log_odds >= 0, "Positive", "Negative")
  
  # Add significance markers
  fisher_results$significance <- ""
  fisher_results$significance[fisher_results$p_adjusted < 0.05] <- "*"
  fisher_results$significance[fisher_results$p_adjusted < 0.01] <- "**"
  fisher_results$significance[fisher_results$p_adjusted < 0.001] <- "***"
  
  return(fisher_results)
}

# Perform Fisher's tests
fisher_results <- perform_fisher_tests(defense_clean %>% filter(type %in% top_defense_systems), 
                                       anti_defense_clean)

# Create matrix for heatmap
# Convert to wide format for the heatmap
heatmap_data <- fisher_results %>%
  select(defense_system, antidefense_system, log_odds, significance) %>%
  # Ensure consistent ordering
  mutate(
    defense_system = factor(defense_system, levels = top_defense_systems),
    antidefense_system = factor(antidefense_system, levels = antidefense_systems)
  )

# Create improved heatmap showing both direction and significance
panel_d <- ggplot(heatmap_data, aes(x = antidefense_system, y = defense_system)) +
  # Use log2 odds ratio for fill color (direction)
  geom_tile(aes(fill = log_odds), color = "white", size = 0.5) +
  
  # Use a diverging color scale centered at zero
  scale_fill_gradient2(
    low = "#2166AC",    # Blue for negative associations
    mid = "white",      # White for neutral
    high = "#B2182B",   # Red for positive associations
    midpoint = 0,
    limits = c(-10, 10),
    name = expression(log[2](OR))
  ) +
  
  # Add significance markers
  geom_text(
    aes(label = significance),
    size = 5
  ) +
  
  # Theming
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 8),
    axis.text.y = element_text(face = "bold", size = 9),
    panel.grid = element_blank(),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10)
  ) +
  
  # Labels
  labs(
    title = "D. Defense-Antidefense Associations",
    subtitle = "* p < 0.05, ** p < 0.01, *** p < 0.001 (FDR-corrected)",
    x = "Anti-Defense System",
    y = "Defense System"
  )

# ==============================================================================
# Combine panels and save results
# ==============================================================================

# Save individual panels
ggsave(file.path(output_dir, "antidefense_panel_a.png"), panel_a, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "antidefense_panel_b.png"), panel_b, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "antidefense_panel_c.png"), panel_c, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "antidefense_panel_d.png"), panel_d, width = 8, height = 6, dpi = 300)

# Combine panels into one figure
combined_plot <- grid.arrange(
  panel_a, panel_b,
  panel_c, panel_d,
  ncol = 2,
  nrow = 2
)

# Save combined figure
ggsave(file.path(output_dir, "Figure10_antidefense_analysis.png"), combined_plot, width = 12, height = 10, dpi = 300)
ggsave(file.path(output_dir, "Figure10_antidefense_analysis.pdf"), combined_plot, width = 12, height = 10)

# Save analysis results
write_csv(fisher_results, file.path(output_dir, "defense_antidefense_correlation_results.csv"))

message("Analysis complete. Results saved to ", output_dir)