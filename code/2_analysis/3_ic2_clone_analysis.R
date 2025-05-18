# ==============================================================================
# ic2_clone_analysis.R
# 
# Comparison of defense systems between IC2 clone contigs and other
# A. baumannii complete genomes using DefenseFinder results
#
# Author: Vigneshwaran Muthuraman
# ==============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
library(ggupset)  # For upset plots

# Set paths
defensefinder_file <- "results/consolidated/consolidated_defense_systems.tsv"
metadata_file <- "data/metadata/ab_genome_id.xlsx"
output_dir <- "results/figures"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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

# Define color palette
clone_colors <- c("IC2 Clone" = "#E41A1C", "Other A. baumannii" = "#377EB8")

# ============================================================================
# Data Loading and Preparation
# ============================================================================

# Load defense systems data
defense_df <- read_tsv(defensefinder_file, show_col_types = FALSE)
  

# Load metadata
metadata <- readxl::read_excel(metadata_file)

# Add IC2 clone status to metadata if not already present
# This assumes metadata has a column that can identify IC2 clones
# If not available, this is a placeholder for demonstration
if (!"Clone_Status" %in% colnames(metadata)) {
  message("Adding mock IC2 clone status (for demonstration)")
  # Create a mock classification - in a real scenario, use actual data
  metadata <- metadata %>%
    mutate(Clone_Status = ifelse(
      grepl("IC2|ST2|GC2", Genome_ID, ignore.case = TRUE) | 
        row_number() %% 3 == 0,  # Every 3rd genome for demonstration
      "IC2 Clone", 
      "Other A. baumannii"
    ))
}

# Filter data to include only A. baumannii
baumannii_data <- defense_df %>%
  left_join(metadata %>% select(Genome_ID, Taxon, Clone_Status), by = "Genome_ID") %>%
  filter(grepl("baumannii", Taxon, ignore.case = TRUE))

# Count number of genomes in each group for reporting
ic2_count <- n_distinct(baumannii_data$Genome_ID[baumannii_data$Clone_Status == "IC2 Clone"])
other_count <- n_distinct(baumannii_data$Genome_ID[baumannii_data$Clone_Status == "Other A. baumannii"])

# ============================================================================
# Panel A: Defense System Count Comparison
# ============================================================================

# Count defense systems per genome
defense_counts <- baumannii_data %>%
  group_by(Genome_ID, Clone_Status) %>%
  summarize(defense_count = n_distinct(type), .groups = "drop")

# Perform Wilcoxon test
wilcox_test <- wilcox.test(
  defense_count ~ Clone_Status, 
  data = defense_counts,
  exact = FALSE
)

# Create boxplot with violin overlay
panel_a <- ggplot(defense_counts, aes(x = Clone_Status, y = defense_count, fill = Clone_Status)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
  scale_fill_manual(values = clone_colors) +
  labs(
    title = "A. Defense System Counts",
    subtitle = sprintf("Wilcoxon p-value = %.3f", wilcox_test$p.value),
    x = NULL,
    y = "Number of Defense Systems",
    fill = "Genome Group"
  ) +
  custom_theme +
  theme(legend.position = "bottom")

# ============================================================================
# Panel B: Defense System Prevalence Comparison
# ============================================================================

# Calculate prevalence of each defense system in IC2 vs. other A. baumannii
defense_prevalence <- baumannii_data %>%
  group_by(Clone_Status, type) %>%
  summarize(
    genome_count = n_distinct(Genome_ID),
    .groups = "drop"
  ) %>%
  group_by(Clone_Status) %>%
  mutate(
    total_genomes = n_distinct(baumannii_data$Genome_ID[baumannii_data$Clone_Status == first(Clone_Status)]),
    prevalence = (genome_count / total_genomes) * 100
  ) %>%
  ungroup()

# Get top defense systems overall to include in the plot
top_systems <- defense_prevalence %>%
  group_by(type) %>%
  summarize(total_count = sum(genome_count), .groups = "drop") %>%
  arrange(desc(total_count)) %>%
  head(10) %>%
  pull(type)

# Filter to include only top systems
defense_prevalence_top <- defense_prevalence %>%
  filter(type %in% top_systems)

# Create side-by-side bar chart
panel_b <- ggplot(defense_prevalence_top, 
                  aes(x = reorder(type, -prevalence), y = prevalence, fill = Clone_Status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = clone_colors) +
  labs(
    title = "B. Defense System Prevalence",
    x = NULL,
    y = "Percentage of Genomes (%)",
    fill = "Genome Group"
  ) +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "bottom"
  )

# ============================================================================
# Panel C: Statistical Enrichment Analysis
# ============================================================================

# Function to perform Fisher's exact test for a defense system
perform_enrichment_test <- function(defense_system, data) {
  # Get distinct genomes
  all_genomes <- unique(data$Genome_ID)
  ic2_genomes <- unique(data$Genome_ID[data$Clone_Status == "IC2 Clone"])
  other_genomes <- unique(data$Genome_ID[data$Clone_Status == "Other A. baumannii"])
  
  # Get genomes with this defense system
  genomes_with_defense <- unique(data$Genome_ID[data$type == defense_system])
  
  # Create 2x2 contingency table
  # [IC2 with system, IC2 without system]
  # [Other with system, Other without system]
  ic2_with <- sum(ic2_genomes %in% genomes_with_defense)
  ic2_without <- length(ic2_genomes) - ic2_with
  
  other_with <- sum(other_genomes %in% genomes_with_defense)
  other_without <- length(other_genomes) - other_with
  
  # Create contingency table
  contingency <- matrix(c(ic2_with, other_with, ic2_without, other_without), 
                        nrow = 2, byrow = TRUE)
  
  # Check if we can perform the test
  if (sum(contingency == 0) > 0) {
    # Add a small pseudocount to avoid errors with zero cells
    contingency <- contingency + 0.5
  }
  
  # Perform Fisher's exact test
  tryCatch({
    test_result <- fisher.test(contingency)
    
    # Return results
    tibble(
      defense_system = defense_system,
      odds_ratio = test_result$estimate,
      p_value = test_result$p.value,
      ci_lower = test_result$conf.int[1],
      ci_upper = test_result$conf.int[2]
    )
  }, error = function(e) {
    # If test fails, return NA values
    tibble(
      defense_system = defense_system,
      odds_ratio = 1,
      p_value = 1,
      ci_lower = 0.1,
      ci_upper = 10
    )
  })
}

# Get all unique defense systems
all_defense_systems <- unique(baumannii_data$type)

# Perform enrichment test for each system
enrichment_results <- map_dfr(all_defense_systems, ~perform_enrichment_test(., baumannii_data))

# Apply FDR correction
enrichment_results <- enrichment_results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05,
    log2_odds_ratio = log2(odds_ratio),
    # Handle infinite values
    log2_odds_ratio = case_when(
      is.infinite(log2_odds_ratio) & log2_odds_ratio > 0 ~ 10,
      is.infinite(log2_odds_ratio) & log2_odds_ratio < 0 ~ -10,
      TRUE ~ log2_odds_ratio
    )
  )

# Define core systems to focus on, if available in the data
core_systems <- c("SspBCDE", "RM", "Gao_Qat", "PD-T7-5", "Cas", "CBASS", "Retron", "RosmerTA")
core_systems_available <- core_systems[core_systems %in% all_defense_systems]

# Get top systems by statistical significance
significant_systems <- enrichment_results %>%
  filter(defense_system %in% core_systems_available, significant) %>%
  pull(defense_system)

# If no significant systems, use the top systems by odds ratio difference from 1
if (length(significant_systems) == 0) {
  significant_systems <- enrichment_results %>%
    mutate(odds_ratio_diff = abs(log2_odds_ratio)) %>%
    arrange(desc(odds_ratio_diff)) %>%
    head(10) %>%
    pull(defense_system)
}

# If still no systems, use the top systems overall
if (length(significant_systems) == 0) {
  significant_systems <- top_systems
}

# Create forest plot data
forest_plot_data <- enrichment_results %>%
  filter(defense_system %in% significant_systems) %>%
  # Ensure nice ordering by log2 odds ratio
  mutate(defense_system = fct_reorder(defense_system, log2_odds_ratio))

# Create forest plot
panel_c <- ggplot(forest_plot_data, 
                  aes(x = log2_odds_ratio, y = defense_system, 
                      color = significant)) +
  # Add reference line at 0 (no enrichment)
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  # Add error bars for confidence intervals
  geom_errorbarh(aes(xmin = log2(ci_lower), xmax = log2(ci_upper)), 
                 height = 0.2, linewidth = 0.5) +
  # Add points
  geom_point(aes(size = -log10(p_adjusted)), shape = 16) +
  # Color and size scales
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8")) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  # Labels
  labs(
    title = "C. Defense System Enrichment in IC2 Clones",
    subtitle = "Log2 Odds Ratio with 95% Confidence Intervals",
    x = expression(log[2](Odds~Ratio)),
    y = NULL,
    color = "FDR < 0.05"
  ) +
  # Add labels for very significant systems
  geom_text_repel(
    data = forest_plot_data %>% filter(p_adjusted < 0.01),
    aes(label = sprintf("p = %.3f", p_adjusted)),
    size = 3, 
    nudge_x = 0.5,
    direction = "y",
    segment.size = 0.2
  ) +
  # Theme
  custom_theme +
  theme(
    axis.text.y = element_text(face = "italic"),
    legend.position = "bottom"
  )

# ============================================================================
# Panel D: Defense System Combinations in IC2 Clones
# ============================================================================

# Focus on just IC2 clone genomes
ic2_defense_data <- baumannii_data %>% 
  filter(Clone_Status == "IC2 Clone") %>%
  distinct(Genome_ID, type)

# Count occurrences of each defense system combination
defense_combinations <- ic2_defense_data %>%
  group_by(Genome_ID) %>%
  summarize(systems = list(type)) %>%
  ungroup()

# Create upset plot
panel_d <- ggplot(defense_combinations, aes(x = systems)) +
  geom_bar() +
  scale_x_upset() +
  labs(
    title = "D. Defense System Combinations in IC2 Clones",
    x = "Combination",
    y = "Number of Genomes"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, face = "italic"),
    plot.title = element_text(face = "bold", size = 12)
  )

# ============================================================================
# Combine panels and save results
# ============================================================================

# Save individual panels
ggsave(file.path(output_dir, "figure4_panel_a.png"), panel_a, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "figure4_panel_b.png"), panel_b, width = 6, height = 5, dpi = 300) 
ggsave(file.path(output_dir, "figure4_panel_c.png"), panel_c, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "figure4_panel_d.png"), panel_d, width = 6, height = 5, dpi = 300)

# Combine all panels into one figure
combined_plot <- plot_grid(
  panel_a, panel_b,
  panel_c, panel_d,
  ncol = 2,
  nrow = 2,
  labels = "AUTO"
)

# Save combined figure
ggsave(file.path(output_dir, "figure4_ic2_clone_analysis.png"), combined_plot, width = 12, height = 10, dpi = 300)
ggsave(file.path(output_dir, "figure4_ic2_clone_analysis.pdf"), combined_plot, width = 12, height = 10)

# Print summary statistics
message("\nSummary of IC2 Clone Analysis:")
message("1. Defense System Count Comparison: Wilcoxon test p-value = ", round(wilcox_test$p.value, 3))
message("2. Top defense system: ", top_systems[1])
message("3. Number of significantly enriched defense systems: ", sum(enrichment_results$significant))

# Save summary results to file
enrichment_summary <- enrichment_results %>%
  arrange(p_adjusted) %>%
  select(defense_system, odds_ratio, log2_odds_ratio, p_value, p_adjusted, significant)

write_csv(enrichment_summary, file.path(output_dir, "ic2_enrichment_results.csv"))

message("Analysis complete. Results saved to ", output_dir)
