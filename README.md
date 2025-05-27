# Defence System Architecture in *Acinetobacter* Species

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Project Overview

This repository contains the code, analysis, and findings from my Master's thesis research on bacterial defence systems in *Acinetobacter* species, with a focus on their relationship to antimicrobial resistance and mobile genetic elements.

*Acinetobacter baumannii* has emerged as a significant healthcare challenge due to its remarkable ability to acquire antimicrobial resistance. This project investigates the complex landscape of bacterial defence systems across different *Acinetobacter* species, examining how these systems might influence the acquisition and maintenance of resistance determinants.

## Key Research Objectives

- Map the distribution and organization of defence systems across *Acinetobacter* species
- Compare defence system predictions between DefenseFinder and PADLOC bioinformatic tools
- Analyze species-specific variations in defence strategies between *A. baumannii*, *A. pittii*, and other *Acinetobacter* species
- Explore correlations between defence systems and antibiotic resistance genes
- Investigate relationships between defence systems and mobile genetic elements

## Dataset

The analysis was performed on two distinct datasets:
- 132 complete *Acinetobacter* genomes (including 43 *A. baumannii*, 27 *A. pittii*, and 62 other *Acinetobacter* species)
- 90 contig-level assemblies of *A. baumannii* clinical isolates


## Repository Structure

- `code/`: R and Bash scripts for computational analysis
- `data/`: Metadata about genomes and example datasets
- `results/`: Key figures and findings
- `docs/`: Additional documentation including thesis summary

## Methodology Overview

The bioinformatic workflow involved:

1. **Genome Analysis**: Using DefenseFinder and PADLOC to predict defence systems
2. **Antibiotic Resistance Gene Detection**: Identification using ResFinder
3. **Mobile Genetic Element Analysis**: BLAST-based identification of integrative mobile elements
4. **Statistical Analysis**: Correlation analysis using Fisher's exact tests and Spearman's rank correlation
5. **Visualization**: Generation of heatmaps, network visualizations, and statistical plots

![Methodology Workflow](results/figures/methodology_workflow.png)

## Prerequisites

- R (version 4.4.0 or higher)
- [DefenseFinder](https://github.com/mdmparis/defense-finder)
- [PADLOC](https://github.com/padlocbio/padloc)
- [ResFinder](https://github.com/genomicepidemiology/resfinder)
- NCBI E-utilities (for genome download)

## Installation

1. Clone this repository: https://github.com/vikos77/Acinetobacter-defence-systems.git

2. Install R and required packages

3. Install bioinformatics tools:
DefenseFinder: Follow instructions at [DefenseFinder](https://github.com/mdmparis/defense-finder)
PADLOC: Follow instructions at [PADLOC](https://github.com/padlocbio/padloc)
Resfinder : Follow instructions at [Resfinder](https://github.com/genomicepidemiology/resfinder)


### Key Findings

### 1. Defense System Distribution and Tool Comparison
![image](results/figures/defense_system_distribution_comparison.png)

**DefenseFinder vs. PADLOC Detection Patterns**:
- DefenseFinder identifies 3-6 defense systems per genome (mean = 5.0), with most genomes clustering around 4 systems
- PADLOC detects higher counts (mean = 8.6), with a broader distribution extending up to 20 systems per genome
- Both tools converge on RM systems as the most prevalent defense mechanism, but show quantitative differences in detection sensitivity

### 2. Restriction-Modification Systems Dominate the Defense Landscape
![image](results/figures/defense_system_prevalence_comparison.png)

**Universal RM Prevalence**:
- RM systems are detected in 102/132 genomes (77.3%) by DefenseFinder, making them the overwhelmingly dominant defense mechanism
- PADLOC identifies RM_type_I in 74/132 genomes (56.1%), confirming RM system prevalence across tools
- RM systems serve as the conserved baseline defense across all *Acinetobacter* species

### 3. Species-Specific Defense System Profiles
![image](results/figures/defense_systems_by_species_simple.png)

**A. baumannii Defense Signature**:

- Enriched in Gao_Qat (17 genomes), SspBCDE (16 genomes), and PD-T7-5 (16 genomes) systems
- SspBCDE emerges as a cornerstone system, frequently occurring as the sole defense mechanism in clinical isolates (47.8% of contig assemblies)

**A. pittii Defense Profile**:
- Characterized by CBASS (10 genomes) and Septu (9 genomes) systems
- Notably lacks SspBCDE enrichment, distinguishing it from *A. baumannii*
- Shows different secondary defense strategies compared to its close relative

**Other Acinetobacter Species**:
- Maintain high RM system prevalence (57 genomes) as baseline defense
- Display more diverse secondary defense profiles including CBASS, RosmerTA, and Cas systems

### 4. Complex Defense System Co-occurrence Networks

#### Network Topology and System Interconnectivity
![image](results/figures/defensefinder_cooccur_circos.png)

**Circos Network Visualization Reveals Defense System Architecture**:
- **RM System Dominance**: The largest sector in the network represents RM systems, with extensive connections to multiple other defense systems, confirming their role as the central hub in the *Acinetobacter* defense network
- **Secondary System Clusters**: Gao_Qat, CBASS, Gabija, and PD-T7-5 form distinct network nodes with varying degrees of connectivity
- **Interconnection Density**: The network shows dense interconnections between certain system pairs, indicating frequent co-occurrence patterns that exceed random expectation
- **Sector Size Proportional to Prevalence**: Visual sector sizes directly correspond to system prevalence (RM: 102 genomes, RosmerTA: 37 genomes, Cas: 30 genomes), providing immediate insight into defense system hierarchy

#### Statistical Co-occurrence Matrix Analysis
![image](results/figures/defensefinder_cooccur_matrix.png)

**Significant Positive Associations (Cooperative Defense Systems)**:
- **Gao_Qat Network**: Shows multiple significant positive associations across the matrix (red circles with yellow backgrounds)
  - Strong co-occurrence with PD-T7-5, PD-T4-5, and multiple secondary systems
  - Suggests functional compatibility and potential synergistic protection mechanisms
- **CBASS-Retron Partnership**: Statistically significant positive association (log₂OR: 2.02, p<0.05) indicates these systems frequently co-occur, possibly targeting different stages of phage infection
- **SspBCDE Associations**: Multiple positive correlations with various defense systems, supporting its role as a cornerstone system that integrates well with other defense mechanisms

**Significant Negative Associations (Mutually Exclusive Systems)**:
- **RM System Exclusivity**: Strong negative correlations (blue circles) with multiple secondary systems:
  - RM vs. PD-T7-5: log₂OR = -3.01 (p<0.05) - suggests functional redundancy or metabolic incompatibility
  - RM vs. Gao_Qat: log₂OR = -3.74 (p<0.05) - indicates alternative evolutionary strategies
- **Matrix Diagonal Significance**: Self-associations (diagonal elements) show maximum significance, validating the analytical approach
- **RloC Antagonism**: Negative correlations with multiple systems suggest this nuclease-based defense may conflict with other mechanisms


### 5. Defense System-Antibiotic Resistance Relationships

![image](results/figures/defence_arg_correlation.png)

#### Hierarchical Clustering Reveals Functional Defense System Groups

**ARG-Facilitative Defense Systems (Top Cluster)**:
The correlation heatmap with hierarchical clustering reveals a distinct cluster of defense systems that consistently show positive associations with antibiotic resistance genes:

- **Gao_Qat Systems**: The most ARG-associated defense system showing significant positive correlations with:
  - Tetracycline resistance: *tetB* (***p<0.001) - strongest association
  - Aminoglycoside resistance: *armA* (***p<0.001), *aph3-Ia* (*p<0.05), *aadA1* (**p<0.01)
  - β-lactam resistance: *blaOXA-66* (*p<0.05), *blaOXA-23* (**p<0.01)
  - Sulfonamide resistance: *sul1* (**p<0.01), *sul2* (*p<0.05)
  - β-lactamase: *blaTDC-25* (*p<0.05)

- **SspBCDE Systems**: Show broad positive associations across multiple ARG classes:
  - Strong correlations with *tetB*, *armA*, *aadA1*, and β-lactamase genes
  - Consistent pattern of ARG facilitation across resistance categories
  - Statistical significance (***p<0.001 to *p<0.05) across multiple associations

- **PD-T7-5 and PD-T4-5 Systems**: Cluster together with moderate positive correlations:
  - Selective associations with specific ARGs rather than broad spectrum
  - Particularly strong with certain β-lactamase and aminoglycoside resistance genes

#### ARG-Restrictive Defense Systems (Bottom Cluster)

**RM Systems: The Primary ARG Barrier**:
- **Contrasting Pattern**: RM systems show the opposite correlation pattern (blue coloring) compared to the ARG-facilitative cluster
- **Broad Negative Associations**: Significant negative correlations with the same ARGs that positively associate with Gao_Qat/SspBCDE
- **Statistical Significance**: Multiple highly significant negative associations (**p<0.01 to ***p<0.001)
- **Functional Interpretation**: Act as barriers to horizontal acquisition of resistance determinants

**Intermediate Defense Systems**:
- **Neutral/Variable Associations**: Systems like Retron, AbiH, Paris show weaker or variable correlations
- **Selective Patterns**: Some systems show positive associations with specific ARG types while being neutral to others
- **Functional Diversity**: Suggests different mechanisms of interaction with mobile genetic elements carrying ARGs

#### Statistical Robustness and Clustering Validation

**Multiple Testing Correction**:
- **FDR Correction Applied**: All p-values adjusted for multiple comparisons using Benjamini-Hochberg procedure
- **Effect Sizes**: Correlation coefficients range from strong negative (-0.7) to strong positive (+0.8)


### 6. Defense System-Mobile Genetic Element Interactions

![image](results/figures/defence_protein_correlation.png)

#### Hierarchical Clustering Reveals Mobile Element Processing Strategies

**IME-Facilitative Defense Systems (Middle Blue Cluster)**:
The correlation analysis reveals a distinct cluster of defense systems that show positive associations with integrative mobile element (IME) proteins, indicating facilitated horizontal gene transfer:

- **SspBCDE Systems: The Primary IME Facilitator**:
  - **Strongest Associations**: Shows highly significant positive correlations (***p<0.001) across all major IME protein functions
  - **Transposase Compatibility**: Strong associations with Tnp25, Tnp26, Tnp28 transposases
  - **Broad Spectrum**: Positive correlations with putative transposases (202-putative, 200-putative) and specialized functions (198-IS26, 196-transposase)
  - **Biological Significance**: Indicates SspBCDE systems actively facilitate or tolerate mobile element integration

- **Gao_Qat Systems: Selective IME Processing**:
  - **Moderate Positive Associations**: Consistent (**p<0.01) positive correlations across IME protein categories
  - **Functional Selectivity**: Particularly strong associations with specific transposase types
  - **Integration with ARG Data**: Matches the pattern seen with antibiotic resistance genes, confirming coordinated facilitation

- **PD-T7-5 and PD-T4-5 Systems: Auxiliary IME Facilitators**:
  - **Coordinated Pattern**: Both systems show similar (*p<0.05) positive association profiles
  - **Functional Clustering**: Hierarchical clustering groups these systems together, indicating similar mechanisms
  - **Complementary Function**: May provide additional layers of mobile element tolerance

#### IME-Restrictive Defense Systems (Top Red and Bottom Pink Clusters)

**Lamassu-Fam and Mokosh: Specialized IME Restrictors**:
- **Strong Negative Associations**: Intense red coloring with high statistical significance (**p<0.01) across most IME protein functions
- **Broad Spectrum Restriction**: Negative correlations with multiple transposase types and mobile element proteins
- **Functional Specialization**: These systems may have evolved specifically to counter mobile genetic element activity
- **Contrasting Strategy**: Opposite approach to SspBCDE cluster, representing alternative evolutionary solutions

**RM Systems: Universal Mobile Element Barriers**:
- **Consistent Negative Pattern**: Significant negative associations (**p<0.01) across all IME protein categories
- **Non-specific Restriction**: Broad-spectrum negative correlations indicate general foreign DNA restriction
- **Classical Function**: Aligns with traditional understanding of RM systems as barriers to horizontal gene transfer
- **Statistical Robustness**: High significance levels across multiple IME protein types

**CBASS Systems: Moderate IME Restriction**:
- **Selective Negative Associations**: Weaker but significant (*p<0.05) negative correlations
- **Intermediate Strategy**: Less restrictive than RM systems but still generally antagonistic to mobile elements

#### IME Protein Functional Categories and Defense System Specificity

**Transposase Enzyme Families**:
- **Tnp25/Tnp26/Tnp28 Series**: Major transposase families showing coordinated association patterns
  - SspBCDE: Strong positive associations across all three
  - RM systems: Consistent negative associations
  - Functional clustering indicates these transposases have similar integration mechanisms

**Putative Transposases (202-putative, 200-putative)**:
- **Uncharacterized Functions**: Defense system associations may help predict functional roles
- **Evolutionary Relationships**: Similar association patterns to characterized transposases suggest related functions
- **Discovery Potential**: Defense system correlations could guide functional annotation efforts

**Specialized Mobile Element Proteins**:
- **IS26 Transposases**: Insertion sequence elements with specific integration preferences
- **196-transposase and 169-transposase**: Distinct transposase families with unique association profiles
- **Functional Diversity**: Different transposase types show varying compatibility with defense systems

### 7. Integrated Genomic Architecture: Cross-Element Correlations

![image](results/figures/final_correlation.png)

#### System-Level Correlation Matrix Reveals Fundamental Genomic Trade-offs

**Defense Systems as Universal Barriers**:
The correlation matrix between genomic element counts reveals defense systems' role as broad-spectrum barriers to foreign genetic material:

- **Defense vs. ARG Relationship**: Weak negative correlation (r = -0.11)
  - Indicates defense systems generally restrict antibiotic resistance gene acquisition
  - Consistent with previous findings showing RM systems' restrictive effects
  - Represents the cost of maintaining comprehensive phage defense

- **Defense vs. IME Relationship**: Weak negative correlation (r = -0.07)
  - Defense systems moderately restrict integrative mobile element accumulation
  - Aligns with detailed IME-defense system analysis showing restriction patterns
  - Suggests defense systems limit horizontal gene transfer vectors

- **Defense vs. Anti-Defense Relationship**: Essentially neutral (r = -0.01)
  - Near-zero correlation indicates independent accumulation patterns
  - May reflect different selective pressures acting on these system types
  - Suggests anti-defense systems operate through different mechanisms than direct interference

#### Mobile Genetic Elements Drive Resistance Dissemination

**ARG-IME Co-mobilization**: Strongest positive correlation (r = 0.69, ***p<0.001)
- **Biological Significance**: Very strong correlation confirms IMEs as primary vehicles for ARG dissemination
- **Clinical Relevance**: Highest correlation coefficient in the matrix indicates tight functional coupling
- **Mechanistic Basis**: Physical linkage of resistance genes on mobile genetic elements
- **Evolutionary Implication**: Co-selection maintains both resistance genes and their mobilization machinery

**Statistical Robustness**:
- **Highest significance level** (***p<0.001) indicates extremely reliable association
- **Large effect size** (r = 0.69) suggests substantial biological importance
- **Clinical prediction**: IME abundance can predict ARG load in clinical isolates

#### Anti-Defense Systems Facilitate Genetic Plasticity

**Anti-Defense as HGT Facilitators**:
- **Anti-Defense vs. ARG**: Moderate positive correlation (r = 0.50, ***p<0.001)
  - Strong statistical significance indicates reliable association
  - Anti-defense systems promote resistance gene acquisition
  - Mechanistic explanation: counteract host defense systems that would otherwise restrict ARG uptake

- **Anti-Defense vs. IME**: Weak-moderate positive correlation (r = 0.26, **p<0.01)
  - Statistically significant association supports facilitatory role
  - Anti-defense systems enable mobile element integration
  - Lower correlation than ARG association suggests indirect relationship

#### Quantitative Framework for Genomic Element Interactions

**Correlation Strength Interpretation**:
- **Strong Positive (r > 0.6)**: ARG-IME relationship represents tight functional coupling
- **Moderate Positive (r = 0.3-0.6)**: Anti-Defense-ARG relationship indicates significant biological interaction
- **Weak Correlations (r < 0.3)**: Defense system relationships suggest subtle but consistent effects
- **Near-Zero Correlations**: Defense-Anti-Defense relationship indicates independent dynamics

**Statistical Significance Hierarchy**:
- **Highest Significance (***p<0.001)**: ARG-IME and Anti-Defense-ARG correlations
- **High Significance (**p<0.01)**: Anti-Defense-IME correlation
- **Non-significant relationships**: All defense system correlations (indicating general restrictive effects rather than specific associations)

#### Evolutionary and Clinical Implications

**Genomic Architecture Models**:

1. **Conservative Genome Architecture**:
   - High defense system counts
   - Low ARG, IME, and anti-defense system counts
   - Emphasis on phage protection over genetic plasticity
   - Represented by RM-dominant genomes

2. **Plastic Genome Architecture**:
   - Moderate defense system counts (selective systems like SspBCDE/Gao_Qat)
   - High ARG and IME counts
   - Elevated anti-defense system presence
   - Optimized for rapid adaptation in clinical environments

**Clinical Prediction Framework**:
- **ARG Load Prediction**: IME count (r = 0.69) is the strongest predictor of resistance gene abundance
- **Genetic Plasticity Assessment**: Anti-defense system count predicts capacity for rapid genetic change
- **Defense Strategy Classification**: Defense system composition indicates evolutionary strategy

#### Mechanistic Integration Across Analyses

**Validation of Previous Findings**:
- **Defense System Dichotomy**: Confirmed by overall negative correlations with foreign genetic elements
- **Mobile Element Centrality**: ARG-IME correlation validates mobile elements as resistance dissemination vehicles
- **Anti-Defense Function**: Positive correlations with both ARGs and IMEs confirm facilitatory role

**Network-Level Understanding**:
- **Hub Elements**: IMEs serve as central hubs connecting resistance genes and facilitating mechanisms
- **Barrier Elements**: Defense systems serve as network barriers reducing overall connectivity
- **Facilitator Elements**: Anti-defense systems increase network connectivity and gene flow

#### Statistical Framework and Methodology

**Spearman Rank Correlation Justification**:
- **Non-parametric approach**: Appropriate for count data without distributional assumptions
- **Robust to outliers**: Resistant to extreme values in genomic element counts
- **Ordinal relationships**: Captures monotonic relationships regardless of linearity

**Multiple Testing Considerations**:
- **Bonferroni correction applied**: Conservative approach for multiple correlation testing

### 8. IC2 Clone-Specific Defense Architecture

**Reduced Defense Complexity**:
- IC2 clones carry significantly fewer defense systems (1-2 per genome) compared to other *A. baumannii* strains (~5 systems)
- SspBCDE system predominates in IC2 clones, often as the sole defense mechanism
- Correlates with enhanced antibiotic resistance gene carriage (average 17.1 ARGs per genome)

**Global Clinical Success**:
- Minimal restriction barriers may facilitate rapid resistance acquisition through horizontal gene transfer
- Trade-off between comprehensive phage defense and genetic plasticity may explain IC2's pandemic success

## Citation

If you use any of the code or findings from this project, please cite:

```
Muthuraman, V. (2025). Defence System Architecture in Acinetobacter: 
Species-Specific Patterns and Implications for Antimicrobial Resistance. 
Masters Thesis, Teesside University.
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgements

I would like to express my sincere gratitude to my supervisor, Dr. Shehreen Saadlee, for her invaluable guidance and support throughout this research project.
