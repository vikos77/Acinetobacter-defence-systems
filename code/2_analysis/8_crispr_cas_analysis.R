# ==============================================================================
# 8_crispr_cas_analysis.R
# 
# Analysis of CRISPR-Cas systems in Acinetobacter species and their
# relationship with other defense systems
# 
# This script uses the existing CRISPR-Cas summary file that contains
# pre-processed CRISPR array and Cas gene counts for each genome
# 
# Author: Vigneshwaran Muthuraman
# ==============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(patchwork)
library(RColorBrewer)

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
crispr_cas_summary_file <- "crispr_cas/CRISPR-Cas_summary.tsv"  # The summary file you showed
metadata_file <- "data/metadata/Acinetobacter_metadata.xlsx"
output_dir <- "crispr_cas/figures"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# Data Loading and Preparation
# ==============================================================================

# Load DefenseFinder data
message("Loading DefenseFinder results...")
defense_df <- read_tsv(defensefinder_file, show_col_types = FALSE)

# Clean DefenseFinder data - ensure unique defense systems per genome
defense_clean <- defense_df %>%
  distinct(Genome_ID, type) %>% 
  mutate(Genome_ID = str_remove(Genome_ID, "\\.[0-9]+$"))

# Load CRISPR-Cas summary data
message("Loading CRISPR-Cas summary data...")
crispr_cas_data <- read_tsv(crispr_cas_summary_file, show_col_types = FALSE)

# Clean column names to remove any special characters
crispr_cas_data <- crispr_cas_data %>%
  rename(
    Genome_ID = `Sequence(s)`,
    CRISPR_arrays = `CRISPR array(s)`,
    Nb_CRISPRs = `Nb CRISPRs`,
    Evidence_levels = `Evidence-levels`,
    Cas_clusters = `Cas cluster(s)`,
    Nb_Cas = `Nb Cas`,
    Cas_Types_Subtypes = `Cas Types/Subtypes`
  )

# Get total genomes from CRISPR-Cas data
all_genomes <- crispr_cas_data$Genome_ID
total_genomes <- length(all_genomes)

message(paste("Total genomes in dataset:", total_genomes))

# ==============================================================================
# CRISPR-Cas System Categorization
# ==============================================================================

# Create comprehensive categorization using the summary data
genome_crispr_status <- crispr_cas_data %>%
  mutate(
    # Convert to numeric, handling NA values
    crispr_count = as.numeric(ifelse(is.na(Nb_CRISPRs) | Nb_CRISPRs == "N/A", 0, Nb_CRISPRs)),
    cas_count = as.numeric(ifelse(is.na(Nb_Cas) | Nb_Cas == "N/A", 0, Nb_Cas)),
    
    # Determine presence/absence
    has_crispr_array = crispr_count > 0,
    has_cas_genes = cas_count > 0,
    
    # Create categories
    crispr_cas_category = case_when(
      has_crispr_array & has_cas_genes ~ "Both CRISPR & Cas",
      has_crispr_array & !has_cas_genes ~ "CRISPR array only",
      !has_crispr_array & has_cas_genes ~ "Cas genes only",
      !has_crispr_array & !has_cas_genes ~ "Neither CRISPR nor Cas"
    )
  ) %>%
  select(Genome_ID, crispr_count, cas_count, has_crispr_array, has_cas_genes, crispr_cas_category)

# Print summary statistics
message("\n=== CRISPR-Cas System Distribution ===")
crispr_summary_stats <- genome_crispr_status %>%
  count(crispr_cas_category, sort = TRUE) %>%
  mutate(percentage = round(n/sum(n) * 100, 1))

print(crispr_summary_stats)

# Show detailed counts
message(paste("Genomes with CRISPR arrays:", sum(genome_crispr_status$has_crispr_array)))
message(paste("Genomes with Cas genes:", sum(genome_crispr_status$has_cas_genes)))
message(paste("Genomes with both:", sum(genome_crispr_status$has_crispr_array & genome_crispr_status$has_cas_genes)))
message(paste("Genomes with neither:", sum(!genome_crispr_status$has_crispr_array & !genome_crispr_status$has_cas_genes)))

# ==============================================================================
# Panel A: CRISPR-Cas System Distribution
# ==============================================================================

# Create distribution bar chart
panel_a <- ggplot(crispr_summary_stats, 
                  aes(x = reorder(crispr_cas_category, n), y = n, fill = crispr_cas_category)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0(n, " (", percentage, "%)")), 
            hjust = -0.1, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "A. CRISPR-Cas System Distribution",
    x = "Category",
    y = "Number of Genomes",
    fill = "Category"
  ) +
  custom_theme +
  theme(legend.position = "none")

# ==============================================================================
# Panel B: Defense System Prevalence Comparison
# ==============================================================================

# Define the three groups for comparison using the updated category names
both_crispr_cas <- genome_crispr_status %>%
  filter(crispr_cas_category == "Both CRISPR & Cas") %>%
  pull(Genome_ID)

neither_crispr_cas <- genome_crispr_status %>%
  filter(crispr_cas_category == "Neither CRISPR nor Cas") %>%
  pull(Genome_ID)

crispr_array_only <- genome_crispr_status %>% 
  filter(crispr_cas_category == "CRISPR array only") %>% 
  pull(Genome_ID)

message("Group sizes:")
message(paste("Both CRISPR & Cas:", length(both_crispr_cas), "genomes"))
message(paste("Neither CRISPR nor Cas:", length(neither_crispr_cas), "genomes"))
message(paste("CRISPR array only:", length(crispr_array_only), "genomes"))

# Calculate defense system prevalence in all three groups
defense_comparison <- defense_clean %>%
  # Include all three groups in the filter
  filter(Genome_ID %in% c(both_crispr_cas, neither_crispr_cas, crispr_array_only)) %>%
  mutate(
    group = case_when(
      Genome_ID %in% both_crispr_cas ~ "Both CRISPR & Cas",
      Genome_ID %in% neither_crispr_cas ~ "Neither CRISPR nor Cas",
      Genome_ID %in% crispr_array_only ~ "CRISPR array only",
      TRUE ~ "Unknown"  # Safety net
    )
  ) %>%
  group_by(group, type) %>%
  summarise(count = n(), .groups = "drop") %>%
  # Add group totals for percentage calculation
  group_by(group) %>%
  mutate(
    total_genomes = case_when(
      group == "Both CRISPR & Cas" ~ length(both_crispr_cas),
      group == "Neither CRISPR nor Cas" ~ length(neither_crispr_cas),
      group == "CRISPR array only" ~ length(crispr_array_only),
      TRUE ~ 1  # Safety net
    ),
    prevalence_percent = round(count / total_genomes * 100, 1)
  ) %>%
  ungroup()

# Get top defense systems (present in at least 3 genomes in any group)
top_defense_systems <- defense_comparison %>%
  group_by(type) %>%
  summarise(max_count = max(count), .groups = "drop") %>%
  filter(max_count >= 3) %>%
  arrange(desc(max_count)) %>%
  head(20) %>%  # Top 20 systems
  pull(type)

# Filter for top systems
defense_comparison_filtered <- defense_comparison %>%
  filter(type %in% top_defense_systems)

# Create grouped bar chart
panel_b <- ggplot(defense_comparison_filtered, 
                  aes(x = reorder(type, prevalence_percent), y = prevalence_percent, fill = group)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(aes(label = paste0(prevalence_percent, "%")), 
            position = position_dodge(width = 0.9), 
            hjust = -0.1, size = 3, family = "Times New Roman") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "B. Defence System Prevalence: CRISPR-Cas Status Comparison",
    subtitle = paste("Comparison between", length(both_crispr_cas), "(Both),", 
                     length(crispr_array_only), "(Array only), and", 
                     length(neither_crispr_cas), "(Neither) genomes"),
    x = "Defence System Type",
    y = "Prevalence (%)",
    fill = "CRISPR-Cas Status"
  ) +
  # Updated color scheme for three groups
  scale_fill_manual(values = c(
    "Both CRISPR & Cas" = "#2E86C1",      # Blue
    "CRISPR array only" = "#F39C12",       # Orange  
    "Neither CRISPR nor Cas" = "#E74C3C"   # Red
  )) +
  theme(
    # All text elements set to Times New Roman
    text = element_text(family = "Times New Roman"),
    axis.text.x = element_text(size = 10, family = "Times New Roman"),
    axis.text.y = element_text(size = 10, family = "Times New Roman"),
    axis.title.x = element_text(size = 12, family = "Times New Roman"),
    axis.title.y = element_text(size = 12, family = "Times New Roman"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14, family = "Times New Roman"),
    plot.subtitle = element_text(size = 11, family = "Times New Roman"),
    legend.title = element_text(size = 10, family = "Times New Roman"),
    legend.text = element_text(size = 9, family = "Times New Roman")
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 5))


print(panel_b)

# ==============================================================================
# Panel C: Statistical Comparison of Defense System Counts
# ==============================================================================

# Count defense systems per genome
defense_counts_per_genome <- defense_clean %>%
  group_by(Genome_ID) %>%
  summarize(defense_count = n(), .groups = "drop")

# Merge with CRISPR-Cas status (ensure consistent genome IDs)
defense_counts_comparison <- genome_crispr_status %>%
  left_join(defense_counts_per_genome, by = "Genome_ID") %>%
  mutate(defense_count = ifelse(is.na(defense_count), 0, defense_count)) %>%
  # Filter for the three main groups
  filter(crispr_cas_category %in% c("Both CRISPR & Cas", "CRISPR array only", "Neither CRISPR nor Cas"))

# Summary statistics
summary_stats <- defense_counts_comparison %>%
  group_by(crispr_cas_category) %>%
  summarise(
    n_genomes = n(),
    mean_defense = round(mean(defense_count), 2),
    median_defense = median(defense_count),
    sd_defense = round(sd(defense_count), 2),
    mean_crispr = round(mean(crispr_count), 2),
    mean_cas = round(mean(cas_count), 2),
    .groups = "drop"
  )

message("\n=== Defense System Count Summary ===")
print(summary_stats)

# Statistical tests
# Kruskal-Wallis test (non-parametric)
kruskal_test <- kruskal.test(defense_count ~ crispr_cas_category, data = defense_counts_comparison)
message(paste("\nKruskal-Wallis test p-value:", format.pval(kruskal_test$p.value, digits = 4)))

# Post-hoc pairwise comparisons if significant
pairwise_comparisons <- NULL
if(kruskal_test$p.value < 0.05) {
  message("\nPerforming post-hoc pairwise comparisons...")
  pairwise_comparisons <- pairwise.wilcox.test(
    defense_counts_comparison$defense_count, 
    defense_counts_comparison$crispr_cas_category, 
    p.adjust.method = "bonferroni"
  )
  print(pairwise_comparisons$p.value)
}

# Create box plot with statistical annotations
panel_c <- ggplot(defense_counts_comparison, 
                  aes(x = crispr_cas_category, y = defense_count, fill = crispr_cas_category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
               fill = "yellow", color = "black") +
  scale_fill_manual(values = c(
    "Both CRISPR & Cas" = "#2E86C1",      # Blue
    "CRISPR array only" = "#F39C12",       # Orange
    "Neither CRISPR nor Cas" = "#E74C3C"   # Red
  )) +
  labs(
    title = "C. Defense System Counts: Box Plot Comparison",
    subtitle = paste("Kruskal-Wallis test p =", format.pval(kruskal_test$p.value, digits = 4)),
    x = "CRISPR-Cas Status", 
    y = "Number of Defense Systems per Genome",
    caption = "Yellow diamonds indicate group means"
  ) +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.caption = element_text(size = 9, color = "gray60")
  )

# Add statistical significance annotations if needed
if (!is.null(pairwise_comparisons)) {
  p_values <- pairwise_comparisons$p.value
  
  # Add significance brackets (customize based on your specific p-values)
  panel_c <- panel_c +
    # Example annotations - adjust coordinates and labels based on your data
    annotate("segment", x = 1, xend = 2, y = 12, yend = 12, color = "black") +
    annotate("text", x = 1.5, y = 12.5, label = "**", size = 5) +
    annotate("segment", x = 1, xend = 3, y = 13.5, yend = 13.5, color = "black") +
    annotate("text", x = 2, y = 14, label = "*", size = 5) +
    annotate("text", x = 2.5, y = 11, label = "ns", size = 3, color = "gray60")
}

# ==============================================================================
# Combine panels and save results
# ==============================================================================

# Create three-panel figure
library(patchwork)

# Define layout manually
layout_design <- "
AB
CB
"

combined_figure <- panel_a + panel_b + panel_c +
  plot_layout(design = layout_design) +
  plot_annotation(
    title = "CRISPR-Cas System Analysis in Acinetobacter genomes",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )
# Save individual panels
ggsave(file.path(output_dir, "crispr_panel_a.png"), panel_a, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "crispr_panel_b.png"), panel_b, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "crispr_panel_c.png"), panel_c, width = 8, height = 6, dpi = 300)

# Save combined figure
ggsave(file.path(output_dir, "Supplement_Figure3_crispr_cas_analysis.png"), 
       combined_figure, width = 16, height = 10, dpi = 300)
ggsave(file.path(output_dir, "Supplement_Figure3_crispr_cas_analysis.pdf"), 
       combined_figure, width = 14, height = 10)

# ==============================================================================
# Additional Analysis: Defense Systems in Genomes Without CRISPR-Cas
# ==============================================================================

# Identify genomes without CRISPR-Cas systems
genomes_without_crispr_cas <- genome_crispr_status %>%
  filter(crispr_cas_category == "Neither CRISPR nor Cas") %>%
  pull(Genome_ID)

message(paste("\n=== Defense Systems in Genomes Without CRISPR-Cas ==="))
message(paste("Analyzing", length(genomes_without_crispr_cas), "genomes without CRISPR-Cas systems..."))

# Analyze other defense systems in these genomes
other_defense_systems <- defense_clean %>%
  filter(Genome_ID %in% genomes_without_crispr_cas) %>%
  count(type, sort = TRUE) %>%
  mutate(percentage = round(n / length(genomes_without_crispr_cas) * 100, 1))

message("Top defense systems in genomes WITHOUT CRISPR-Cas:")
print(head(other_defense_systems, 10))

# ==============================================================================
# Save Analysis Results
# ==============================================================================

# Save summary data
write_csv(genome_crispr_status, file.path(output_dir, "crispr_cas_genome_classification.csv"))
write_csv(crispr_summary_stats, file.path(output_dir, "crispr_cas_summary_statistics.csv"))
write_csv(defense_prevalence_filtered, file.path(output_dir, "defense_prevalence_by_crispr_status.csv"))
write_csv(other_defense_systems, file.path(output_dir, "defense_systems_no_crispr_cas.csv"))
write_csv(summary_stats, file.path(output_dir, "defense_count_summary_by_crispr.csv"))

# Save statistical test results
if (!is.null(pairwise_comparisons)) {
  # Convert pairwise comparison results to a data frame
  pairwise_df <- data.frame(
    comparison = rownames(pairwise_comparisons$p.value),
    p_value = as.vector(pairwise_comparisons$p.value),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(p_value))
  
  write_csv(pairwise_df, file.path(output_dir, "crispr_defense_pairwise_comparisons.csv"))
}

message("\n=== Analysis Complete ===")
message(paste("Results saved to:", output_dir))
message("Key findings:")
message(paste("1. Genomes with both CRISPR array and Cas genes:", 
              sum(genome_crispr_status$crispr_cas_category == "Both CRISPR & Cas")))
message(paste("2. Genomes with CRISPR array only:", 
              sum(genome_crispr_status$crispr_cas_category == "CRISPR array only")))
message(paste("3. Genomes with neither CRISPR nor Cas:", 
              sum(genome_crispr_status$crispr_cas_category == "Neither CRISPR nor Cas")))
message(paste("4. Statistical significance (Kruskal-Wallis p-value):", 
              format.pval(kruskal_test$p.value, digits = 4)))