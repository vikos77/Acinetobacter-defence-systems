# Detailed Methodology

This document provides a comprehensive overview of the bioinformatic and statistical methods used in this project.

## Bioinformatic Workflow

### 1. Data Acquisition and Preparation

Two primary datasets were used in this analysis:

**Dataset 1**: 100 *Acinetobacter baumannii* genomes
- 90 contig-level assemblies
- 10 complete genome assemblies

**Dataset 2**: 122 complete *Acinetobacter* spp. chromosomes
- Downloaded from NCBI using accession numbers

For the primary analysis, we used a combined dataset of 132 complete genomes (10 from Dataset 1 and all 122 from Dataset 2).

### 2. Computational Environment

- **Operating System**: Ubuntu 24.04.1 LTS
- **Package Management**: Conda (miniconda) version 24.11.3
- **R Environment**: R version 4.4.2 and RStudio version 4.4.2

### 3. Defence System Prediction

Two specialized bioinformatic tools were used:

#### DefenseFinder v1.2.0

```bash
defense-finder run -o results/defensefinder_output/[GenomeID] [GenomeFASTA]
```

Key output files:
- `[GenomeID]_defence_finder_genes.tsv`: Gene-level data
- `[GenomeID]_defence_finder_hmmer.tsv`: HMMER hit data
- `[GenomeID]_defence_finder_systems.tsv`: System-level summary

#### PADLOC v1.0.1

Two-step process:
1. CRISPR array detection with CRISPRCasFinder
2. PADLOC analysis with CRISPR integration

```bash
padloc --genome [GenomeFASTA] --out results/padloc_output/[GenomeID] --crisp [CRISPROutputGFF]
```

### 4. Anti-Defence System Prediction

```bash
defense-finder run --antidefensefinder-only -o results/antidefense_output/[GenomeID] [GenomeFASTA]
```

### 5. Data Consolidation

Custom bash scripts were used to consolidate outputs from both tools into tab-separated value files for downstream analysis.

## Statistical Analysis

### 1. Data Processing and Deduplication

To ensure accurate analysis, duplicate entries for each defence system type within the same genome were removed:

```R
defense_df_deduplicated <- defense_df %>%
  distinct(Genome_ID, type, .keep_all = TRUE)
```

### 2. Correlation Analyses

Multiple correlation methods were employed:

- **Spearman's Rank Correlation**: Used for examining relationships between defence system counts, anti-defence system counts, IME counts, and ARG counts
- **Fisher's Exact Test**: For assessing associations between specific defence system pairs, defence-ARG pairs, and defence-IME pairs
- **Multiple Testing Correction**: Benjamini-Hochberg procedure for controlling false discovery rate
- **Co-occurrence Analysis**: Jaccard similarity indices for quantifying co-occurrence patterns

### 3. Visualization Techniques

Various visualization methods were employed:

#### Distribution Analysis
```R
ggplot(defense_counts, aes(x = count)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Defence Systems per Genome",
       x = "Number of Defence Systems",
       y = "Number of Genomes")
```

#### Correlation Heatmaps
```R
pheatmap(correlation_matrix,
         color = colorRampPalette(c("navy", "white", "firebrick"))(100),
         breaks = seq(-1, 1, length.out = 101),
         display_numbers = display_matrix,
         main = "Correlation Between Defense Systems and ARGs")
```

#### Network Visualization
For defense system co-occurrence network visualization, we utilized custom network analysis functions.

## Analysis of Defence-Resistance Relationships

To investigate relationships between defence systems and antibiotic resistance genes (ARGs):

1. **Binary Matrix Creation**: Created presence/absence matrices for both defence systems and ARGs
2. **Fisher's Exact Test**: Applied to each defence system-ARG pair
3. **Odds Ratio Calculation**: Log2-transformed odds ratios were calculated to quantify association strength
4. **Visualization**: Results were visualized using bar charts and heatmaps

## Reproducibility

All analyses were conducted using documented scripts with fixed random seeds to ensure reproducibility. The entire workflow can be executed using the main pipeline script:

```bash
bash code/pipeline/main_pipeline.sh
```

This will reproduce all analyses from raw genome data to final visualizations.
