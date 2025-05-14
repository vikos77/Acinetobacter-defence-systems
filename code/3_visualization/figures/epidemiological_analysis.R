#=================================================================
# Epidemiological_analysis.R
# 
# Creates a multi-panel figure showing epidemiological analysis of
# Acinetobacter dataset:
#   Panel A: Geographic distribution of isolates
#   Panel B: Isolation source distribution
#   Panel C: Genome size vs GC content correlation
#   Panel D: Genomic feature comparison across species groups
#=================================================================
#Load the required libraries
library(ggplot2)
library(maps)
library(patchwork) 
library(ggrepel)
library(countrycode)
library(scales)


# Read metadata
metadata <- readxl::read_xlsx("Acinetobacter_metadata.xlsx")


# Clean up data
metadata <- metadata %>%
  mutate(
    GenomeSize_Mb = `Size(Bp)` / 1000000,
    Species = case_when(
      grepl("baumannii", Taxon) ~ "A. baumannii",
      grepl("pittii", Taxon) ~ "A. pittii",
      grepl("calcoaceticus", Taxon) ~ "A. calcoaceticus",
      grepl("nosocomialis", Taxon) ~ "A. nosocomialis",
      grepl("haemolyticus", Taxon) ~ "A. haemolyticus",
      grepl("junii", Taxon) ~ "A. junii",
      grepl("lwoffi", Taxon) ~ "A. lwoffii",
      TRUE ~ "Other Acinetobacter spp."
    ),
    Country = case_when(
      grepl(":", Location) ~ trimws(str_split_fixed(Location, ":", 2)[,1]),
      !is.na(Location) ~ Location,
      TRUE ~ "Unknown"
    ),
    SourceType = case_when(
      grepl("blood|sputum|urine|swab|aspirate|abscess|wound|Hip|Mouth|Homo|homo", Source, ignore.case=TRUE) ~ "Clinical - Human",
      grepl("feces|stool|Dung", Source, ignore.case=TRUE) ~ "Fecal",
      grepl("soil|water|sewage|waste|recycle|hospital|wastewater|mud|permafrost", Source, ignore.case=TRUE) ~ "Environmental",
      grepl("animal|fish|chicken|cow|pig|Horse|goose", Source, ignore.case=TRUE) ~ "Animal-associated",
      TRUE ~ "Other/Unknown"
    )
  )


# Count isolates by country for map
country_counts <- metadata %>%
  filter(!is.na(Country) & Country != "Unknown") %>%
  count(Country, name = "IsolateCount") %>%
  arrange(desc(IsolateCount))

# Get top countries for labeling
top_countries <- country_counts %>% head(6)

top_countries <- data.frame(
  Country = c("China", "Mexico", "USA", "Germany", "Australia", "Chile"),
  longitude = c(104.2, -102.5, -95.7, 10.4, 133.8, -71.0),
  latitude = c(35.9, 23.6, 37.1, 51.2, -25.3, -35.7),
  IsolateCount = c(39, 13, 12, 7, 5, 4)
)


# Get world map data
world_map <- map_data("world")


map_plot_simple <- ggplot() +
  # Base world map
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white", size = 0.1) +
  # Add points for countries with isolates
  geom_point(data = top_countries,
             aes(x = longitude, y = latitude, size = IsolateCount),
             color = "red", alpha = 0.7) +
  # Add labels for top countries
  geom_label_repel(data = top_countries,
                   aes(x = longitude, y = latitude, 
                       label = paste0(Country, " (", IsolateCount, ")")),
                   size = 3, box.padding = 0.5, point.padding = 0.5,
                   segment.color = "black") +
  scale_size_continuous(range = c(3, 10), name = "Number of\nIsolates") +
  theme_minimal() +
  labs(title = "A. Geographic Distribution") +
  theme(axis.title = element_blank(),
        panel.grid = element_blank())


# Panel B: Genome characteristics

genome_plot <- ggplot(metadata, aes(x = `Size(Bp)`, y = `GC content (%)`)) +
  # Add density contours
  geom_density_2d(color = "gray70") +
  # Add points with uniform color
  geom_point(color = "#3366CC", alpha = 0.7) +
  # Add overall trend line
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(hjust = 1, face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5)
  )+
  scale_x_continuous(
    name   = "Genome size (Mb)",
    labels = function(x) x / 1000000
  ) +
  labs(title = "C. Genome Characteristics", 
       x = "Genome Size (Kb)", 
       y = "GC Content (%)")

# Panel C: Source distribution

source_bar <- metadata %>%
  count(SourceType) %>%
  mutate(
    prop = n / sum(n),
    percentage = round(prop * 100),
    label = paste0(percentage, "%"),
    # Ensure consistent ordering
    SourceType = factor(SourceType, 
                        levels = c("Other/Unknown", "Clinical - Human", 
                                   "Environmental", "Fecal", "Animal-associated"))
  ) %>%
  ggplot(aes(x = reorder(SourceType, n), y = n, fill = SourceType)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = label), hjust = -0.3, size = 5, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  # Set appropriate x-axis limit to accommodate labels
  scale_y_continuous(limits = c(0, max(metadata %>% count(SourceType) %>% pull(n)) * 1.15)) +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.position = "none", 
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 14),
    plot.title = element_text(size = 12)
  ) +
  labs(title = "B. Isolation Sources", 
       x = "", 
       y = "Number of Genomes")

#Panel D

# Compare key metrics across three species groups: A. baumannii, A. pittii, and Others
species_comparison <- metadata %>%
  mutate(SpeciesGroup = case_when(
    Species == "A. baumannii" ~ "Acinetobacter baumannii",
    Species == "A. pittii" ~ "Acinetobacter pittii",
    TRUE ~ "Other Acinetobacter spp."
  )) %>%
  group_by(SpeciesGroup) %>%
  summarize(
    `Mean Genome Size (Kb)` = round(mean(`Size(Bp)`/1000, na.rm = TRUE), 2),
    `Mean GC Content (%)` = round(mean(`GC content (%)`, na.rm = TRUE), 1),
    `Mean Protein Coding Genes` = round(mean(`RefSeq(Protein coding)`, na.rm = TRUE), 0),
    `Genomes Analyzed` = n()
  ) %>%
  # Ensure consistent ordering
  mutate(SpeciesGroup = factor(SpeciesGroup, 
                               levels = c("Acinetobacter baumannii", 
                                          "Acinetobacter pittii", 
                                          "Other Acinetobacter spp."))) %>%
  pivot_longer(cols = -SpeciesGroup, names_to = "Metric", values_to = "Value") %>%
  # Ensure consistent metric ordering
  mutate(Metric = factor(Metric, 
                         levels = c("Genomes Analyzed", 
                                    "Mean Genome Size (Kb)", 
                                    "Mean GC Content (%)", 
                                    "Mean Protein Coding Genes",
                                    "Coding Density (genes/Mb)")))

# Create alternative display labels for the plot
species_labels <- c(
  "Acinetobacter baumannii" = "A. baumannii",
  "Acinetobacter pittii" = "A. pittii",
  "Other Acinetobacter spp." = "Other spp."
)

# Create the heatmap with custom species name labels
species_heatmap <- ggplot(species_comparison, 
                          aes(x = SpeciesGroup, y = Metric, fill = Value)) +
  geom_tile(color = "white", size = 0.8) +
  geom_text(aes(label = Value), size = 3.8, fontface = "bold") +
  scale_fill_gradient(low = "#E6F0FF", high = "#3366CC") +
  # Use custom labels that was created for the x-axis
  scale_x_discrete(labels = species_labels) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "italic", size = 14),
    axis.text.y = element_text(hjust = 1, face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5)
  ) +
  labs(title = "D. Genomic Feature Comparison",
       x = "Species Group") +
  theme(legend.position = "none")


# Combine all plots
combined_plot <- (map_plot_simple + source_bar) / (genome_plot + species_heatmap) +
  plot_layout(guides = "collect") & theme(legend.position = "none")

# Save the plot
ggsave("figures/Figure1_Epidemiological_Analysis.pdf", combined_plot, width = 12, height = 10, dpi = 300)
ggsave("figures/Figure1_Epidemiological_Analysis.png", combined_plot, width = 12, height = 10, dpi = 300)

#Individualpanels
ggsave("figures/Figure1_panel_A.png", map_plot_simple, width = 12, height = 10, dpi = 300)
ggsave("figures/Figure1_panel_B.png", genome_plot, width = 12, height = 10, dpi = 300)
ggsave("figures/Figure1_panel_C.png", source_bar, width = 12, height = 12, dpi = 300)
ggsave("figures/Figure1_panel_D.png", species_heatmap, width = 12, height = 10, dpi = 300)
  
