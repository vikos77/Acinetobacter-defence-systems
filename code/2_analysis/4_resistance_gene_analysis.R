# ==============================================================================
# 4_resistance_gene_analysis.R
# 
# Analysis of antibiotic resistance genes in Acinetobacter species and
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
resfinder_file <- "results/consolidated/consolidated_resfinder_results.tsv"
defensefinder_file <- "results/consolidated/consolidated_defense_systems.tsv"
output_dir <- "results/figures"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# Load and prepare data
# ==============================================================================

# Load ResFinder data
resfinder_df <- read_tsv(resfinder_file, show_col_types = FALSE)

# Load DefenseFinder data for correlation analysis
defense_df <- read_tsv(defensefinder_file, show_col_types = FALSE)

# Clean ResFinder data - focus on unique resistance gene occurrences
resfinder_clean <- resfinder_df %>%
  mutate(Genome_ID = str_extract(Contig, "^([^ ]+)")) %>%  # Extract text before the first space as Genome_ID
  # Assuming ResFinder file has columns like Genome_ID, Resistance gene, etc.
  # Standardize resistance gene column name
  rename_with(~ ifelse(. == "Resistance gene", "Resistance_gene", .), 
              matches("Resistance gene|Resistance_gene")) %>%
  select(Genome_ID, Resistance_gene) %>%
  # Ensure unique resistance gene per genome
  distinct(Genome_ID, Resistance_gene) 

# Clean DefenseFinder data - ensure unique defense systems per genome
defense_clean <- defense_df %>%
  mutate(Genome_ID = gsub("\\.(fas|fasta)$", "", Genome_ID)) %>% 
  distinct(Genome_ID, type) 

# Count total genomes in dataset
total_genomes <- n_distinct(c(resfinder_clean$Genome_ID, defense_clean$Genome_ID))
message(paste("Total unique genomes in dataset:", total_genomes))

# ==============================================================================
# Panel A: Distribution of ARG Counts per Genome
# ==============================================================================

# Count ARGs per genome
args_per_genome <- resfinder_clean %>%
  count(Genome_ID, name = "arg_count") %>%
  # Add genomes with zero ARGs, if any
  right_join(data.frame(Genome_ID = unique(c(resfinder_clean$Genome_ID, defense_clean$Genome_ID))), 
             by = "Genome_ID") %>%
  mutate(arg_count = ifelse(is.na(arg_count), 0, arg_count))

# Create histogram of ARG counts
panel_a <- ggplot(args_per_genome, aes(x = arg_count)) +
  geom_histogram(binwidth = 1, fill = "#E41A1C", color = "black", alpha = 0.8) +
  geom_vline(aes(xintercept = mean(arg_count, na.rm = TRUE)), 
             linetype = "dashed", color = "darkred", size = 1) +
  annotate("text", x = mean(args_per_genome$arg_count) + 1, 
           y = max(table(args_per_genome$arg_count)) * 0.9, 
           label = paste("Mean:", round(mean(args_per_genome$arg_count), 1)),
           color = "darkred", fontface = "bold") +
  labs(
    title = "A. Distribution of Antibiotic Resistance Gene Counts",
    x = "Number of ARGs per Genome",
    y = "Number of Genomes"
  ) +
  custom_theme

# ==============================================================================
# Panel B: Prevalence of Different Resistance Genes
# ==============================================================================

# Count and calculate prevalence of each resistance gene
arg_prevalence <- resfinder_clean %>%
  count(Resistance_gene, name = "count") %>%
  arrange(desc(count)) %>%
  mutate(
    percentage = count / total_genomes * 100,
    Resistance_gene = fct_reorder(Resistance_gene, count, .desc = TRUE)
  )

# Get top 20 resistance genes for visualization
top_args <- arg_prevalence %>%
  head(20)

# Create bar chart
panel_b <- ggplot(top_args, aes(x = Resistance_gene, y = count)) +
  geom_bar(stat = "identity", fill = "#E41A1C", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "B. Top 20 Antibiotic Resistance Genes",
    x = NULL,
    y = "Number of Genomes"
  ) +
  custom_theme +
  theme(axis.text.y = element_text(face = "italic", size = 9))

# ==============================================================================
# Panel C: Resistance Gene Categories
# ==============================================================================

# Categorize resistance genes
categorize_args <- function(gene) {
  if (grepl("blaOXA|blaTEM|blaADC|blaNDM|blaVIM|blaIMP|blaCTX", gene)) return("β-lactamases")
  if (grepl("aac|aad|aph|ant|str", gene)) return("Aminoglycosides")
  if (grepl("tet", gene)) return("Tetracyclines")
  if (grepl("sul", gene)) return("Sulfonamides")
  if (grepl("dfr", gene)) return("Trimethoprim")
  if (grepl("mph|erm|msr", gene)) return("Macrolides")
  if (grepl("qnr|gyr|par", gene)) return("Quinolones")
  if (grepl("cat|cml|floR", gene)) return("Chloramphenicol")
  return("Other")
}

# Add category to each resistance gene
arg_with_category <- resfinder_clean %>%
  mutate(Category = sapply(Resistance_gene, categorize_args))

# Summarize by category
arg_categories <- arg_with_category %>%
  count(Category, name = "count") %>%
  arrange(desc(count)) %>%
  mutate(
    percentage = count / sum(count) * 100,
    Category = fct_reorder(Category, count, .desc = TRUE)
  )

# Create category pie chart

panel_c <- ggplot(arg_categories, aes(x = 2, y = count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, position = position_stack(vjust = 0.95)) +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "C. Resistance Gene Categories",
    fill = "Category"
  ) +
  theme_void() +
  xlim(0.5, 2.5) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  ) +
  geom_text(
    aes(label = paste0(round(percentage), "%")),
    position = position_stack(vjust = 0.5),
    size = 2
  )
# ==============================================================================
# Panel D: Correlation Between Defense Systems and ARGs
# ==============================================================================

# Get top defense systems for analysis
top_defense_systems <- defense_clean %>%
  count(type, sort = TRUE) %>%
  head(20) %>%
  pull(type)

# Get top ARGs for analysis
top_args_list <- arg_prevalence %>%
  head(10) %>%
  pull(Resistance_gene)

# Create presence-absence matrices
# Defense systems matrix
defense_matrix <- defense_clean %>%
  filter(type %in% top_defense_systems) %>%
  mutate(present = 1) %>%
  pivot_wider(
    id_cols = Genome_ID,
    names_from = type,
    values_from = present,
    values_fill = 0
  )

# ARG matrix
arg_matrix <- resfinder_clean %>%
  filter(Resistance_gene %in% top_args_list) %>%
  mutate(present = 1) %>%
  pivot_wider(
    id_cols = Genome_ID,
    names_from = Resistance_gene,
    values_from = present,
    values_fill = 0
  )

# Join matrices
all_genomes <- unique(c(defense_matrix$Genome_ID, arg_matrix$Genome_ID))
joined_matrix <- data.frame(Genome_ID = all_genomes) %>%
  left_join(defense_matrix, by = "Genome_ID") %>%
  left_join(arg_matrix, by = "Genome_ID") %>%
  mutate(across(everything(), ~replace_na(., 0)))

# Extract defense systems and ARGs columns
defense_cols <- top_defense_systems[top_defense_systems %in% colnames(joined_matrix)]
arg_cols <- top_args_list[top_args_list %in% colnames(joined_matrix)]

# Initialize result matrices
fisher_results <- data.frame(
  Defense_System = character(),
  ARG = character(),
  OddsRatio = numeric(),
  Pvalue = numeric(),
  stringsAsFactors = FALSE
)

# Perform Fisher's exact test for each defense system - ARG pair
for (defense in defense_cols) {
  for (arg in arg_cols) {
    defense_data <- joined_matrix[[defense]]
    arg_data <- joined_matrix[[arg]]
    
    contingency_table <- table(
      factor(defense_data, levels = c(0, 1)),
      factor(arg_data, levels = c(0, 1))
    )
    
    # Add small constant if table has zeros
    if (any(contingency_table == 0)) {
      contingency_table <- contingency_table + 0.5
    }
    
    # Perform Fisher's exact test
    test_result <- fisher.test(contingency_table)
    
    # Store results
    fisher_results <- rbind(fisher_results, data.frame(
      Defense_System = defense,
      ARG = arg,
      OddsRatio = test_result$estimate,
      Pvalue = test_result$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

# Apply FDR correction
fisher_results <- fisher_results %>%
  mutate(
    Padj = p.adjust(Pvalue, method = "BH"),
    Significant = Padj < 0.05,
    LogOddsRatio = log2(OddsRatio),
    # Handle infinite values
    LogOddsRatio = case_when(
      is.infinite(LogOddsRatio) & LogOddsRatio > 0 ~ 16,
      is.infinite(LogOddsRatio) & LogOddsRatio < 0 ~ -16,
      TRUE ~ LogOddsRatio
    )
  )

# Create matrix for heatmap
heatmap_matrix <- matrix(
  fisher_results$LogOddsRatio,
  nrow = length(defense_cols),
  ncol = length(arg_cols),
  dimnames = list(defense_cols, arg_cols)
)

# Create significance indicator matrix
significance_matrix <- matrix(
  "", 
  nrow = length(defense_cols),
  ncol = length(arg_cols),
  dimnames = list(defense_cols, arg_cols)
)

# Fill significance indicators
for (i in 1:nrow(fisher_results)) {
  if (fisher_results$Padj[i] < 0.05) {
    if (fisher_results$Padj[i] < 0.001) {
      significance_matrix[fisher_results$Defense_System[i], fisher_results$ARG[i]] <- "***"
    } else if (fisher_results$Padj[i] < 0.01) {
      significance_matrix[fisher_results$Defense_System[i], fisher_results$ARG[i]] <- "**"
    } else {
      significance_matrix[fisher_results$Defense_System[i], fisher_results$ARG[i]] <- "*"
    }
  }
}

# Create temp directory if needed
temp_dir <- file.path(output_dir, "temp_arg_analysis")
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE)
}

# Ensure all graphics devices are closed
while (!is.null(dev.list())) {
  dev.off()
}

# LARGER PANEL PARAMETERS
# Increase dimensions for better visibility in multi-panel figure
cellsize <- 40        # Increased from 30
fontsize_row <- 18    # Increased from 16
fontsize_col <- 18    # Increased from 16
fontsize_main <- 22   # Increased from 18
fontsize_num <- 18    # Increased from 16

italic_col_labels <- parse(text = paste0(
  "italic(\"", colnames(heatmap_matrix), "\")"
))

# Create heatmap using pheatmap
pdf(file.path(output_dir, "resistance_defense_correlation_heatmap.pdf"), width = 16, height = 14)
pheatmap(
  heatmap_matrix,
  display_numbers = significance_matrix,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  breaks = seq(-10, 10, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = fontsize_row,
  fontsize_col = fontsize_col,
  cellwidth = cellsize,
  cellheight = cellsize,
  main = "Correlation Between Defence Systems and ARGs",
  fontsize = fontsize_main,
  number_color = "black",
  fontsize_number = fontsize_num,
  border_color = NA,
  angle_col = 45,
  labels_col = italic_col_labels,
)
dev.off()

# Save heatmap as PNG
png(file.path(output_dir, "resistance_defense_correlation_heatmap.png"), width = 16, height = 14, units = "in", res = 300)
pheatmap(
  heatmap_matrix,
  display_numbers = significance_matrix,
  color = colorRampPalette(c("navy", "white", "firebrick"))(100),
  breaks = seq(-10, 10, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = fontsize_row,
  fontsize_col = fontsize_col,
  cellwidth = cellsize,
  cellheight = cellsize,
  main = "Correlation Between Defence Systems and ARGs",
  fontsize = fontsize_main,
  number_color = "black",
  fontsize_number = fontsize_num,
  border_color = NA,
  angle_col = 45,
  labels_col = italic_col_labels,
)
dev.off()

# Save panel D as a ggplot (alternative to pheatmap)
fisher_long <- fisher_results %>%
  mutate(
    Defense_System = factor(Defense_System, levels = rev(defense_cols)),
    ARG = factor(ARG, levels = arg_cols)
  )

panel_d <- ggplot(fisher_long, aes(x = ARG, y = Defense_System, fill = LogOddsRatio)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "navy", 
    mid = "white", 
    high = "firebrick",
    midpoint = 0,
    name = expression(log[2](OR))
  ) +
  geom_text(aes(label = ifelse(Significant, 
                               ifelse(Padj < 0.001, "***", 
                                      ifelse(Padj < 0.01, "**", "*")), 
                               "")),
            color = "black", size = 4) +
  labs(
    title = "D. Correlation Between Defence Systems and ARGs",
    x = "Antibiotic Resistance Genes",
    y = "Defence Systems"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )

# ==============================================================================
# Combine panels and save results
# ==============================================================================

# Save individual panels
ggsave(file.path(output_dir, "arg_panel_a.png"), panel_a, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "arg_panel_b.png"), panel_b, width = 6, height = 6, dpi = 300)
ggsave(file.path(output_dir, "arg_panel_c.png"), panel_c, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "arg_panel_d.png"), panel_d, width = 8, height = 6, dpi = 300)

# Combine panels into one figure
combined_plot <- grid.arrange(
  panel_a, panel_b,
  panel_c, panel_d,
  ncol = 2,
  nrow = 2
)

# Save combined figure
ggsave(file.path(output_dir, "Figure9_resistance_analysis.png"), combined_plot, width = 12, height = 10, dpi = 300)
ggsave(file.path(output_dir, "Figure9_resistance_analysis.pdf"), combined_plot, width = 12, height = 10)

# Save analysis results
write_csv(fisher_results, file.path(output_dir, "defense_arg_correlation_results.csv"))