# ==============================================================================
# HMRG_multipanel_analysis.R
# 
# Comprehensive analysis of Heavy Metal Resistance Genes (HMRGs) in 
# Acinetobacter species - combining Fisher's exact tests, species-specific
# patterns, and isolation source comparisons
#
# This script combines three analysis components:
# 1. Defense system-HMRG correlation analysis (Fisher's exact tests)
# 2. Species-specific ARG-HMRG patterns  
# 3. Clinical vs Environmental source comparison
#
# Author: Vigneshwaran Muthuraman
# ==============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(patchwork)

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
hmrg_blast_file <- "Test/HMRG_proteins_vs_acinetobacter.tblastn"
defensefinder_file <- "Test/defense_df.csv"
resfinder_file <- "Test/arg_data.csv"
metadata_file <- "Test/Acinetobacter_genomes_sem2_2025.xlsx"
output_dir <- "Test/figures"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# Data loading and preparation
# ==============================================================================

# Load HMRG BLAST data
hmrg_blast_data <- read.delim(hmrg_blast_file, header = FALSE)

# Rename columns
col_names <- c("QueryID", "SubjectID", "PercentIdentity", "AlignmentLength", 
               "Mismatches", "GapOpens", "QueryStart", "QueryEnd", 
               "SubjectStart", "SubjectEnd", "Evalue", "BitScore",
               "QueryLength", "QueryCovS", "QueryCovHSP")

names(hmrg_blast_data)[1:min(ncol(hmrg_blast_data), length(col_names))] <- col_names

# Filter hits by quality (manuscript criteria)
filtered_hits <- hmrg_blast_data %>%
  filter(PercentIdentity >= 80, QueryCovS >= 80, Evalue <= 0.005)

# Parse HMRG gene names and genome IDs
parsed_hits <- filtered_hits %>%
  mutate(
    # Extract genome ID from subject sequence ID
    genome_id = str_extract(SubjectID, "^\\S+"),
    
    # Extract HMRG gene name from query ID (BacMet format)
    gene_name = case_when(
      grepl("\\|", QueryID) ~ {
        parts <- str_split_fixed(QueryID, "\\|", 5)
        toupper(parts[, 2])
      },
      TRUE ~ toupper(QueryID)
    )
  ) %>%
  # Remove duplicates: same genome with same HMRG gene
  distinct(genome_id, gene_name, .keep_all = TRUE)

message(paste("Parsed", nrow(parsed_hits), "unique HMRG hits"))
message(paste("Unique genomes with HMRGs:", n_distinct(parsed_hits$genome_id)))
message(paste("Unique HMRG genes:", n_distinct(parsed_hits$gene_name)))

# Load defense systems data
defense_df <- read_csv(defensefinder_file, show_col_types = FALSE) %>%
  distinct(Genome_ID, type)

# Load ARG data
arg_df <- read_csv(resfinder_file, show_col_types = FALSE)

# Load metadata
metadata <- readxl::read_excel(metadata_file)

# ==============================================================================
# Panel A: Defense system-HMRG correlation heatmap
# ==============================================================================

# Get top defense systems and HMRG genes
top_defense_systems <- defense_df %>%
  count(type, sort = TRUE) %>%
  head(20) %>%
  pull(type)

top_hmrg_genes <- parsed_hits %>%
  count(gene_name, sort = TRUE) %>%
  head(10) %>%
  pull(gene_name)

# Create binary matrices
defense_matrix <- defense_df %>%
  filter(type %in% top_defense_systems) %>%
  distinct(Genome_ID, type) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = type, values_from = present, values_fill = 0)

hmrg_matrix <- parsed_hits %>%
  filter(gene_name %in% top_hmrg_genes) %>%
  distinct(genome_id, gene_name) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = gene_name, values_from = present, values_fill = 0)

# Get all unique genomes and merge matrices
all_genomes <- unique(c(defense_matrix$Genome_ID, hmrg_matrix$genome_id))
combined_matrix <- data.frame(Genome_ID = all_genomes) %>%
  left_join(defense_matrix, by = "Genome_ID") %>%
  left_join(hmrg_matrix, by = c("Genome_ID" = "genome_id")) %>%
  mutate(across(everything(), ~replace_na(., 0))) %>%
  mutate(across(-Genome_ID, as.numeric))

# Perform Fisher's exact tests
fisher_results <- data.frame()

for (defense in top_defense_systems) {
  for (hmrg in top_hmrg_genes) {
    if (defense %in% colnames(combined_matrix) && hmrg %in% colnames(combined_matrix)) {
      
      # Create 2x2 contingency table
      n11 <- sum(combined_matrix[[defense]] == 1 & combined_matrix[[hmrg]] == 1)
      n10 <- sum(combined_matrix[[defense]] == 1 & combined_matrix[[hmrg]] == 0)
      n01 <- sum(combined_matrix[[defense]] == 0 & combined_matrix[[hmrg]] == 1)
      n00 <- sum(combined_matrix[[defense]] == 0 & combined_matrix[[hmrg]] == 0)
      
      contingency_table <- matrix(c(n11, n01, n10, n00), nrow = 2, byrow = TRUE)
      
      # Add small constant if zeros present
      if (any(contingency_table == 0)) {
        contingency_table <- contingency_table + 0.5
      }
      
      # Fisher's exact test
      tryCatch({
        test_result <- fisher.test(contingency_table)
        fisher_results <- rbind(fisher_results, data.frame(
          Defense_System = defense,
          HMRG_Gene = hmrg,
          OddsRatio = test_result$estimate,
          Pvalue = test_result$p.value
        ))
      }, error = function(e) {
        message("Error in Fisher's test for ", defense, " vs ", hmrg)
      })
    }
  }
}

# Apply FDR correction and create visualization data
fisher_results <- fisher_results %>%
  mutate(
    Padj = p.adjust(Pvalue, method = "BH"),
    LogOddsRatio = log2(OddsRatio),
    LogOddsRatio = case_when(
      is.infinite(LogOddsRatio) & LogOddsRatio > 0 ~ 10,
      is.infinite(LogOddsRatio) & LogOddsRatio < 0 ~ -10,
      TRUE ~ LogOddsRatio
    ),
    Significance = case_when(
      Padj < 0.001 ~ "***",
      Padj < 0.01 ~ "**", 
      Padj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Create matrices for pheatmap
heatmap_data <- matrix(NA, nrow = length(top_defense_systems), ncol = length(top_hmrg_genes))
rownames(heatmap_data) <- top_defense_systems
colnames(heatmap_data) <- top_hmrg_genes

sig_matrix <- matrix("", nrow = length(top_defense_systems), ncol = length(top_hmrg_genes))
rownames(sig_matrix) <- top_defense_systems
colnames(sig_matrix) <- top_hmrg_genes

# Fill matrices with data
for (i in 1:nrow(fisher_results)) {
  defense <- fisher_results$Defense_System[i]
  hmrg <- fisher_results$HMRG_Gene[i]
  heatmap_data[defense, hmrg] <- fisher_results$LogOddsRatio[i]
  sig_matrix[defense, hmrg] <- fisher_results$Significance[i]
}

# Generate pheatmap heatmap (following original HMRG analysis pattern)
ph <- pheatmap(
  heatmap_data,
  display_numbers = sig_matrix,
  color            = colorRampPalette(c("blue","white","firebrick"))(100),
  breaks           = seq(-10, 10, length.out = 101),
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  fontsize_row     = 16,
  fontsize_col     = 12,
  cellwidth        = 45,
  cellheight       = 40,
  main             = "C. Correlation of defence systems and HMRG",    
  fontsize         = 16,
  number_color     = "black",
  fontsize_number  = 16,
  border_color     = "black",
  angle_col        = 45,
  silent           = TRUE
)

# Ensure any existing graphics devices are closed
while (!is.null(dev.list())) {
  dev.off()
}

# Open PNG device
png(file.path(output_dir, "defense_hmrg_heatmap.png"),
    width  = 10*300, height = 14*300, res = 300)

# Draw the heatmap
grid.newpage()
grid.draw(ph$gtable)

# Add "Defence systems" label on the right
grid.text(
  "Defence systems",
  x    = unit(0.92, "npc"),   # near right edge
  y    = unit(0.5,  "npc"),   # centered vertically
  rot  = -90,                 # vertical orientation
  gp   = gpar(fontsize = 16, fontface = "bold")
)

# Add "Heavy metal resistance genes" label at the bottom
grid.text(
  "Heavy metal resistance genes",
  x    = unit(0.47,  "npc"),   # centered horizontally
  y    = unit(0.01, "npc"),   # near bottom edge
  gp   = gpar(fontsize = 16, fontface = "bold")
)

# Close device
dev.off()

# Check if heatmap file was created successfully
if (!file.exists(file.path(output_dir, "defense_hmrg_heatmap.png"))) {
  stop("Failed to create heatmap PNG file")
}

# Load the generated heatmap PNG as grob for combined figure
panel_a <- grid::rasterGrob(
  png::readPNG(file.path(output_dir, "defense_hmrg_heatmap.png")),
  interpolate = TRUE
)
# ==============================================================================
# Panel B: Species-specific analysis with pairwise comparisons
# ==============================================================================

# Prepare data by species
arg_counts <- arg_df %>%
  select(Genome_ID, `Resistance gene`) %>%
  distinct() %>%
  group_by(Genome_ID) %>%
  summarize(arg_count = n(), .groups = "drop")

hmrg_counts <- parsed_hits %>%
  group_by(genome_id) %>%
  summarize(hmrg_count = n(), .groups = "drop")

species_data <- metadata %>%
  mutate(
    Species_Group = case_when(
      grepl("baumannii", Taxon, ignore.case = TRUE) ~ "A. baumannii",
      grepl("pittii", Taxon, ignore.case = TRUE) ~ "A. pittii",
      TRUE ~ "Other Acinetobacter spp."
    )
  ) %>%
  left_join(arg_counts, by = c("Accession Number" = "Genome_ID")) %>%
  left_join(hmrg_counts, by = c("Accession Number" = "genome_id")) %>%
  mutate(
    arg_count = ifelse(is.na(arg_count), 0, arg_count),
    hmrg_count = ifelse(is.na(hmrg_count), 0, hmrg_count)
  )

# Statistical tests
kruskal_arg <- kruskal.test(arg_count ~ Species_Group, data = species_data)
kruskal_hmrg <- kruskal.test(hmrg_count ~ Species_Group, data = species_data)

# Pairwise comparisons with Bonferroni correction
pairwise_arg <- pairwise.wilcox.test(species_data$arg_count, 
                                     species_data$Species_Group, 
                                     p.adjust.method = "bonferroni")
pairwise_hmrg <- pairwise.wilcox.test(species_data$hmrg_count, 
                                      species_data$Species_Group, 
                                      p.adjust.method = "bonferroni")

# Prepare plot data
plot_data <- species_data %>%
  select(`Accession Number`, Species_Group, arg_count, hmrg_count) %>%
  pivot_longer(cols = c(arg_count, hmrg_count), 
               names_to = "Gene_Type", values_to = "Count") %>%
  mutate(Gene_Type = case_when(
    Gene_Type == "arg_count" ~ "ARGs",
    Gene_Type == "hmrg_count" ~ "HMRGs"
  ))

# Function to get significance stars
get_stars <- function(p_value) {
  if (is.na(p_value)) return("ns")
  if (p_value < 0.001) return("***")
  if (p_value < 0.01) return("**")
  if (p_value < 0.05) return("*")
  return("ns")
}

# Extract p-values for annotations
max_arg <- max(species_data$arg_count, na.rm = TRUE)
max_hmrg <- max(species_data$hmrg_count, na.rm = TRUE)

# Create species comparison plot
panel_b <- ggplot(plot_data, aes(x = Species_Group, y = Count, fill = Species_Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  scale_fill_manual(values = c(
    "A. baumannii" = "#E74C3C",
    "A. pittii" = "#3498DB", 
    "Other Acinetobacter spp." = "#2ECC71"
  )) +
  facet_wrap(~ Gene_Type, scales = "free_y") +
  
  # Add pairwise significance annotations
  geom_segment(data = data.frame(Gene_Type = c("ARGs", "HMRGs"), 
                                 x1 = 1, x2 = 2, 
                                 y = c(max_arg * 1.05, max_hmrg * 1.05)),
               aes(x = x1, xend = x2, y = y, yend = y), 
               inherit.aes = FALSE, color = "black") +
  geom_text(data = data.frame(
    Gene_Type = c("ARGs", "HMRGs"),
    x = 1.5, 
    y = c(max_arg * 1.07, max_hmrg * 1.07),
    label = c(
      get_stars(pairwise_arg$p.value["A. pittii", "A. baumannii"]),
      get_stars(pairwise_hmrg$p.value["A. pittii", "A. baumannii"])
    )),
    aes(x = x, y = y, label = label), 
    inherit.aes = FALSE, size = 3, fontface = "bold") +
  
  geom_segment(data = data.frame(Gene_Type = c("ARGs", "HMRGs"), 
                                 x1 = 1, x2 = 3, 
                                 y = c(max_arg * 1.25, max_hmrg * 1.25)),
               aes(x = x1, xend = x2, y = y, yend = y), 
               inherit.aes = FALSE, color = "black") +
  geom_text(data = data.frame(
    Gene_Type = c("ARGs", "HMRGs"),
    x = 2, 
    y = c(max_arg * 1.27, max_hmrg * 1.27),
    label = c(
      get_stars(pairwise_arg$p.value["Other Acinetobacter spp.", "A. baumannii"]),
      get_stars(pairwise_hmrg$p.value["Other Acinetobacter spp.", "A. baumannii"])
    )),
    aes(x = x, y = y, label = label), 
    inherit.aes = FALSE, size = 3, fontface = "bold") +
  
  geom_segment(data = data.frame(Gene_Type = c("ARGs", "HMRGs"), 
                                 x1 = 2, x2 = 3, 
                                 y = c(max_arg * 1.15, max_hmrg * 1.15)),
               aes(x = x1, xend = x2, y = y, yend = y), 
               inherit.aes = FALSE, color = "black") +
  geom_text(data = data.frame(
    Gene_Type = c("ARGs", "HMRGs"),
    x = 2.5, 
    y = c(max_arg * 1.17, max_hmrg * 1.17),
    label = c(
      get_stars(pairwise_arg$p.value["Other Acinetobacter spp.", "A. pittii"]),
      get_stars(pairwise_hmrg$p.value["Other Acinetobacter spp.", "A. pittii"])
    )),
    aes(x = x, y = y, label = label), 
    inherit.aes = FALSE, size = 3, fontface = "bold") +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  labs(
    title = "A. Species-Specific Resistance Gene Patterns",
    subtitle = "Pairwise comparisons: * p<0.05, ** p<0.01, *** p<0.001, ns = not significant",
    x = "Species Group",
    y = "Gene Count per Genome"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35)))

# ==============================================================================
# Panel C: Clinical vs Environmental comparison
# ==============================================================================

# Create source categories
source_data <- metadata %>%
  mutate(
    SourceType = case_when(
      grepl("blood|sputum|urine|swab|aspirate|abscess|wound|Hip|Mouth|Homo|homo", Source, ignore.case = TRUE) ~ "Clinical",
      grepl("soil|water|sewage|waste|recycle|hospital|wastewater|mud|permafrost", Source, ignore.case = TRUE) ~ "Environmental",
      TRUE ~ "Other"
    )
  ) %>%
  filter(SourceType %in% c("Clinical", "Environmental")) %>%
  left_join(arg_counts, by = c("Accession Number" = "Genome_ID")) %>%
  left_join(hmrg_counts, by = c("Accession Number" = "genome_id")) %>%
  mutate(
    arg_count = ifelse(is.na(arg_count), 0, arg_count),
    hmrg_count = ifelse(is.na(hmrg_count), 0, hmrg_count)
  )

# Statistical tests for clinical vs environmental
if (nrow(source_data) > 0 && length(unique(source_data$SourceType)) == 2) {
  clinical_data <- source_data %>% filter(SourceType == "Clinical")
  env_data <- source_data %>% filter(SourceType == "Environmental")
  
  wilcox_arg <- wilcox.test(clinical_data$arg_count, env_data$arg_count)
  wilcox_hmrg <- wilcox.test(clinical_data$hmrg_count, env_data$hmrg_count)
  
  # Prepare plot data
  source_plot_data <- source_data %>%
    select(`Accession Number`, SourceType, arg_count, hmrg_count) %>%
    pivot_longer(cols = c(arg_count, hmrg_count), 
                 names_to = "Gene_Type", values_to = "Count") %>%
    mutate(Gene_Type = case_when(
      Gene_Type == "arg_count" ~ "ARGs",
      Gene_Type == "hmrg_count" ~ "HMRGs"
    ))
  
  max_arg_src <- max(source_data$arg_count, na.rm = TRUE)
  max_hmrg_src <- max(source_data$hmrg_count, na.rm = TRUE)
  
  # Create source comparison plot
  panel_c <- ggplot(source_plot_data, aes(x = SourceType, y = Count, fill = SourceType)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 1) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    scale_fill_manual(values = c("Clinical" = "#E67E22", "Environmental" = "#27AE60")) +
    facet_wrap(~ Gene_Type, scales = "free_y") +
    
    # Add significance annotations
    geom_segment(data = data.frame(Gene_Type = c("ARGs", "HMRGs"), 
                                   x1 = 1, x2 = 2, 
                                   y = c(max_arg_src * 1.1, max_hmrg_src * 1.1)),
                 aes(x = x1, xend = x2, y = y, yend = y), 
                 inherit.aes = FALSE, color = "black") +
    geom_text(data = data.frame(
      Gene_Type = c("ARGs", "HMRGs"),
      x = 1.5, 
      y = c(max_arg_src * 1.12, max_hmrg_src * 1.12),
      label = c(get_stars(wilcox_arg$p.value), get_stars(wilcox_hmrg$p.value))
    ),
    aes(x = x, y = y, label = label), 
    inherit.aes = FALSE, size = 4, fontface = "bold") +
    
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "none"
    ) +
    labs(
      title = "B. Clinical vs Environmental Comparison",
      x = "Isolation Source",
      y = "Gene Count per Genome"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
} else {
  panel_c <- ggplot() + 
    geom_text(aes(x = 0.5, y = 0.5, label = "No clinical/environmental data available"), 
              size = 5) +
    theme_void() +
    labs(title = "Clinical vs Environmental Comparison")
}

# ==============================================================================
# Combine panels and save results
# ==============================================================================

# Create 2-panel layout using grid.arrange since panel_a is a raster grob

left_panels_grob <- arrangeGrob(panel_b, panel_c, ncol = 1, heights = c(1, 1))

# Create combined figure using grid.arrange
combined_figure <- arrangeGrob(
   left_panels_grob, panel_a,
  ncol = 2, 
  widths = c(1, 1.2),
  top = textGrob("Heavy Metal Resistance Gene (HMRG) Analysis in Acinetobacter", 
                 gp = gpar(fontsize = 16, fontface = "bold"))
)

# Save individual panels
# Copy the heatmap file to the standard panel naming
file.copy(file.path(output_dir, "defense_hmrg_heatmap.png"), 
          file.path(output_dir, "hmrg_panel_a_correlation.png"), 
          overwrite = TRUE)

ggsave(file.path(output_dir, "hmrg_panel_b_species.png"), panel_b, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, "hmrg_panel_c_source.png"), panel_c, width = 8, height = 6, dpi = 300)


# Ensure any existing graphics devices are closed
while (!is.null(dev.list())) {
  dev.off()
}

png(file.path(output_dir, "Supplement_Figure4_HMRG_analysis.png"), 
    width = 18*300, height = 12*300, res = 300)
grid.newpage()
grid.draw(combined_figure)
dev.off()

pdf(file.path(output_dir, "Supplement_Figure4_HMRG_analysis.pdf"), 
    width = 18, height = 12)
grid.newpage()
grid.draw(combined_figure)
dev.off()

# Save analysis results
write_csv(fisher_results, file.path(output_dir, "hmrg_defense_correlation_results.csv"))

# Save species summary
species_summary <- species_data %>%
  group_by(Species_Group) %>%
  summarize(
    n_genomes = n(),
    mean_arg = round(mean(arg_count), 2),
    sd_arg = round(sd(arg_count), 2),
    mean_hmrg = round(mean(hmrg_count), 2),
    sd_hmrg = round(sd(hmrg_count), 2),
    .groups = "drop"
  )

write_csv(species_summary, file.path(output_dir, "hmrg_species_summary.csv"))

# Save pairwise results
if (exists("pairwise_arg")) {
  pairwise_summary <- data.frame(
    Comparison = c(
      "A. baumannii vs A. pittii",
      "A. baumannii vs Other Acinetobacter spp.",
      "A. pittii vs Other Acinetobacter spp."
    ),
    ARG_p_value = c(
      pairwise_arg$p.value["A. pittii", "A. baumannii"],
      pairwise_arg$p.value["Other Acinetobacter spp.", "A. baumannii"],
      pairwise_arg$p.value["Other Acinetobacter spp.", "A. pittii"]
    ),
    HMRG_p_value = c(
      pairwise_hmrg$p.value["A. pittii", "A. baumannii"],
      pairwise_hmrg$p.value["Other Acinetobacter spp.", "A. baumannii"],
      pairwise_hmrg$p.value["Other Acinetobacter spp.", "A. pittii"]
    )
  )
  write_csv(pairwise_summary, file.path(output_dir, "hmrg_species_pairwise_comparisons.csv"))
}

message("=== HMRG Analysis Complete ===")
message("Results saved to: ", output_dir)
message("Main outputs:")
message("- Supplement_Figure4_HMRG_analysis.png (2-panel layout)")
message("- Individual panel files and CSV summaries")
