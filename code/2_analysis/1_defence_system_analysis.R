# ==============================================================================
# defence_system_analysis.R
# 
# Analysis of defence systems in Acinetobacter species using
# consolidated outputs from DefenseFinder and PADLOC
#
# Author: Vigneshwaran Muthuraman
# ==============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
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

#-----------------------------------------------------------------
# Data Loading and Preparation
#-----------------------------------------------------------------

# Set file paths based on the repository structure
defensefinder_file <- "results/consolidated/consolidated_defense_systems.tsv"
padloc_file <- "results/consolidated/consolidated_padloc_results.tsv"
metadata_file <- "data/metadata/genome_metadata.xlsx" 

# Create output directory for figures
output_dir <- "results/figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load DefenseFinder data

defense_df <- read_tsv(defensefinder_file, show_col_types = FALSE)

# Clean DefenseFinder data - ensure unique defence systems per genome
defense_df_clean <- defense_df %>%
  group_by(Genome_ID) %>%
  distinct(type, .keep_all = TRUE) %>%
  ungroup()

# Load PADLOC data

padloc_df <- read_tsv(padloc_file, show_col_types = FALSE)

# Clean PADLOC data - ensure unique defence systems per genome
padloc_df_clean <- padloc_df %>%
  group_by(Genome_ID) %>%
  distinct(system, .keep_all = TRUE) %>%
  ungroup()

# Load metadata with taxonomy information

metadata <- readxl::read_excel(metadata_file)  

# Count total genomes analyzed
total_genomes_defense <- n_distinct(defense_df_clean$Genome_ID)
total_genomes_padloc <- n_distinct(padloc_df_clean$Genome_ID)

#-----------------------------------------------------------------
# Panel A: Distribution of Defence System Counts per Genome
#-----------------------------------------------------------------

# Count systems per genome for DefenseFinder
systems_per_genome_defense <- defense_df_clean %>%
  group_by(Genome_ID) %>%
  summarize(num_systems = n()) %>%
  ungroup()

# Count systems per genome for PADLOC
systems_per_genome_padloc <- padloc_df_clean %>%
  group_by(Genome_ID) %>%
  summarize(num_systems = n()) %>%
  ungroup()

# Create histogram for DefenseFinder
panel_a_defense <- ggplot(systems_per_genome_defense, aes(x = num_systems)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#377EB8", alpha = 0.8) +
  labs(
    title = "DefenseFinder",
    subtitle = paste("Mean:", round(mean(systems_per_genome_defense$num_systems), 1)),
    x = "Number of Defence Systems per Genome",
    y = "Number of Genomes"
  ) +
  scale_x_continuous(breaks = seq(0, max(systems_per_genome_defense$num_systems), by = 1)) +
  custom_theme

# Create histogram for PADLOC
panel_a_padloc <- ggplot(systems_per_genome_padloc, aes(x = num_systems)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#E41A1C", alpha = 0.8) +
  labs(
    title = "PADLOC",
    subtitle = paste("Mean:", round(mean(systems_per_genome_padloc$num_systems), 1)),
    x = "Number of Defence Systems per Genome",
    y = "Number of Genomes"
  ) +
  scale_x_continuous(breaks = seq(0, max(systems_per_genome_padloc$num_systems), by = 1)) +
  custom_theme

# Combine DefenseFinder and PADLOC histograms
panel_a <- grid.arrange(
  panel_a_defense, panel_a_padloc,
  ncol = 2,
  top = textGrob("A. Distribution of Defence System Counts", 
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

#-----------------------------------------------------------------
# Panel B: Prevalence of Defense System Types
#-----------------------------------------------------------------

# Count occurrences of each defence system type for DefenseFinder
defense_type_counts_defense <- defense_df_clean %>%
  count(type, name = "count") %>%
  arrange(desc(count)) %>%
  mutate(percentage = count / total_genomes_defense * 100)

# Count occurrences of each defence system type for PADLOC
defense_type_counts_padloc <- padloc_df_clean %>%
  count(system, name = "count") %>%
  arrange(desc(count)) %>%
  mutate(percentage = count / total_genomes_padloc * 100)

# Create bar chart for DefenseFinder (top 20)
panel_b_defense <- defense_type_counts_defense %>%
  head(20) %>%
  ggplot(aes(x = count, y = reorder(type, count))) +
  geom_bar(stat = "identity", fill = "#377EB8", alpha = 0.8) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", count, percentage)), 
            hjust = -0.1, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "DefenseFinder",
    x = "Number of Genomes",
    y = NULL
  ) +
  custom_theme +
  theme(axis.text.y = element_text(face = "bold", size = 10))

# Create bar chart for PADLOC (top 20)
panel_b_padloc <- defense_type_counts_padloc %>%
  head(20) %>%
  ggplot(aes(x = count, y = reorder(system, count))) +
  geom_bar(stat = "identity", fill = "#E41A1C", alpha = 0.8) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", count, percentage)), 
            hjust = -0.1, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "PADLOC",
    x = "Number of Genomes",
    y = NULL
  ) +
  custom_theme +
  theme(axis.text.y = element_text(face = "bold", size = 10))

# Combine both prevalence plots
panel_b <- grid.arrange(
  panel_b_defense, panel_b_padloc,
  ncol = 2,
  top = textGrob("B. Top 20 Defence System Types", 
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

#-----------------------------------------------------------------
# Panel C: Species-Specific Defence System Comparison
#-----------------------------------------------------------------

# Join defence systems with metadata for species information
defense_with_species <- defense_df_clean %>%
  left_join(metadata %>% select(Genome_ID, Taxon), by = "Genome_ID")

padloc_with_species <- padloc_df_clean %>%
  left_join(metadata %>% select(Genome_ID, Taxon), by = "Genome_ID")

# Create species groups
defense_species_grouped <- defense_with_species %>%
  mutate(Species_Group = case_when(
    grepl("baumannii", Taxon, ignore.case = TRUE) ~ "Acinetobacter baumannii",
    grepl("pittii", Taxon, ignore.case = TRUE) ~ "Acinetobacter pittii",
    TRUE ~ "Other Acinetobacter spp."
  ))

padloc_species_grouped <- padloc_with_species %>%
  mutate(Species_Group = case_when(
    grepl("baumannii", Taxon, ignore.case = TRUE) ~ "Acinetobacter baumannii",
    grepl("pittii", Taxon, ignore.case = TRUE) ~ "Acinetobacter pittii",
    TRUE ~ "Other Acinetobacter spp."
  ))

# Function to create species-specific bar charts
create_species_bar <- function(data, species_name, system_col, tool_name, n_top = 10, color = "blue") {
  # Filter data for specific species and get top N systems
  species_data <- data %>%
    filter(Species_Group == species_name) %>%
    count(!!sym(system_col), name = "count") %>%
    arrange(desc(count)) %>%
    head(n_top)
  
  # Calculate max value for proper x-axis scaling
  max_count <- max(species_data$count, na.rm = TRUE)
  x_limit <- ceiling(max_count * 1.15)
  
  # Create horizontal bar chart
  ggplot(species_data, aes(x = count, y = reorder(!!sym(system_col), count))) +
    geom_bar(stat = "identity", fill = color, alpha = 0.8) +
    geom_text(aes(label = count), hjust = -0.2, size = 3) +
    scale_x_continuous(limits = c(0, x_limit)) +
    labs(
      title = paste(species_name, "-", tool_name),
      x = "Number of Genomes",
      y = NULL
    ) +
    custom_theme +
    theme(
      plot.title = element_text(hjust = 0.5, face = "italic", size = 11),
      axis.text.y = element_text(face = "bold", size = 9)
    )
}

# Create plots for each species and tool combination
# A. baumannii
baumannii_defense <- create_species_bar(defense_species_grouped, "Acinetobacter baumannii", 
                                        "type", "DefenseFinder", color = "#377EB8")
baumannii_padloc <- create_species_bar(padloc_species_grouped, "Acinetobacter baumannii", 
                                       "system", "PADLOC", color = "#E41A1C")

# A. pittii
pittii_defense <- create_species_bar(defense_species_grouped, "Acinetobacter pittii", 
                                     "type", "DefenseFinder", color = "#377EB8")
pittii_padloc <- create_species_bar(padloc_species_grouped, "Acinetobacter pittii", 
                                    "system", "PADLOC", color = "#E41A1C")

# Other species
others_defense <- create_species_bar(defense_species_grouped, "Other Acinetobacter spp.", 
                                     "type", "DefenseFinder", color = "#377EB8")
others_padloc <- create_species_bar(padloc_species_grouped, "Other Acinetobacter spp.", 
                                    "system", "PADLOC", color = "#E41A1C")

# Combine all species plots
panel_c <- grid.arrange(
  baumannii_defense, baumannii_padloc,
  pittii_defense, pittii_padloc,
  others_defense, others_padloc,
  ncol = 2, nrow = 3,
  top = textGrob("C. Top Defence Systems by Species Group", 
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

#-----------------------------------------------------------------
# Save individual panels and combined figure
#-----------------------------------------------------------------

# Save Panel A
ggsave(file.path(output_dir, "figure3_panel_a.png"), panel_a, 
       width = 10, height = 5, dpi = 300)

# Save Panel B
ggsave(file.path(output_dir, "figure3_panel_b.png"), panel_b, 
       width = 12, height = 6, dpi = 300)

# Save Panel C
ggsave(file.path(output_dir, "figure3_panel_c.png"), panel_c, 
       width = 10, height = 12, dpi = 300)

