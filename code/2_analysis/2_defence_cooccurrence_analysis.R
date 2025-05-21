  # ==============================================================================
  # defence_cooccurrence_analysis.R
  # 
  # Co-occurrence analysis of defense systems in Acinetobacter using
  # consolidated outputs from DefenseFinder and PADLOC
  #
  # Author: Vigneshwaran Muthuraman
  # ==============================================================================
  
  # Load required libraries
  library(tidyverse)
  library(patchwork)
  library(circlize)
  library(grid)
  library(gridExtra)
  library(png)
  library(viridis)
  
  # Set paths
  defensefinder_file <- "results/consolidated//consolidated_defense_systems.tsv"
  padloc_file <- "results/consolidated//consolidated_padloc_results.tsv"
  output_dir <- "results/figures"
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create temporary directory for intermediate files
  temp_dir <- file.path(output_dir, "temp")
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }
  
  # Define custom theme for consistent visualization
  custom_theme <- theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )
  
  
  # Load data
  defense_df <- read_tsv(defensefinder_file, show_col_types = FALSE)
  padloc_df <- read_csv(padloc_file, show_col_types = FALSE)
  
  # Clean data - keep only unique defence system types per genome
  defense_df_clean <- defense_df %>%
    group_by(Genome_ID) %>%
    distinct(type, .keep_all = TRUE) %>%
    ungroup()
  
  padloc_df_clean <- padloc_df %>%
    group_by(Genome_ID) %>%
    distinct(system, .keep_all = TRUE) %>%
    ungroup()
  
  padloc_df_simplified <- padloc_df_clean %>%
    mutate(
      # Create a simplified system type by extracting main category
      Simple_System_Type = case_when(
        # Handle RM systems (all types)
        grepl("^RM_", system) ~ "RM",
        # Handle CBASS systems
        grepl("^cbass_", system) ~ "CBASS",
        # Handle cas systems
        grepl("^cas_", system) ~ "Cas",
        # Handle DRT systems
        grepl("^DRT_", system) ~ "DRT",
        # Handle septu systems
        grepl("^septu_", system) ~ "Septu",
        # Handle CRISPR
        grepl("CRISPR", system) ~ "CRISPR",
        #Handle PDC systems
        grepl("^PDC-", system) ~ "PDC",
        # Keep qatABCD as is (similar to Gao_Qat in DefenseFinder)
        system == "qatABCD" ~ "Gao_Qat",
        # Keep other types as is
        TRUE ~ system
      )
    ) %>% 
    group_by(Genome_ID, Simple_System_Type) %>%
    slice(1) %>%  # Keep just one entry per genome and simplified system type
    ungroup()
  
  
  # ========== DEFENSEFINDER CO-OCCURRENCE ANALYSIS ==========
  
  # Function to perform co-occurrence analysis and create visualizations
  analyze_cooccurrence <- function(data, system_col, tool_name, colors_palette, output_prefix) {
    # Count occurrences of each defense system
    system_counts <- data %>%
      count(!!sym(system_col), sort = TRUE) %>%
      mutate(prevalence = n / n_distinct(data$Genome_ID) * 100)
    
    # Select top systems for analysis
    top_systems <- system_counts %>%
      head(20) %>%
      pull(!!sym(system_col))
    
    # Create presence-absence matrix
    presence_matrix <- data %>%
      filter(!!sym(system_col) %in% top_systems) %>%
      distinct(Genome_ID, !!sym(system_col)) %>%
      mutate(
        present = 1,
        !!sym(system_col) := factor(!!sym(system_col), levels = top_systems)
      ) %>%
      pivot_wider(
        names_from = !!sym(system_col),
        values_from = present,
        values_fill = 0
      ) %>%
      column_to_rownames("Genome_ID")
    
    # Dynamically generate colors for top systems
    # Use a visually distinct palette like viridis
    system_colors_dynamic <- viridis::viridis_pal(option = "turbo")(length(top_systems))
    names(system_colors_dynamic) <- top_systems
    system_colors_subset <- system_colors_dynamic

    
    # ========== PART 1: CIRCOS PLOT ==========
    
    # Calculate co-occurrence counts
    cooccur_counts <- crossprod(as.matrix(presence_matrix))
    
    # Apply threshold for better visualization
    threshold <- max(cooccur_counts) * 0.05
    filtered_counts <- cooccur_counts
    filtered_counts[filtered_counts < threshold] <- 0
    
    # Create circos plot
    circos_file <- file.path(temp_dir, paste0(output_prefix, "_circos.png"))
    png(circos_file, width = 2000, height = 1600, res = 200)
    
    # Set up circos plot
    circos.clear()
    circos.par(gap.after = 1)
    par(mar = c(0.5, 0.5, 0.5, 0.5))
    
    # Plot chord diagram
    chordDiagram(
      filtered_counts,
      grid.col = system_colors_subset,
      transparency = 0.5,
      directional = 0,
      annotationTrack = c("grid", "axis"),
      preAllocateTracks = list(track.height = 0.1)
    )
    
    # Add labels
    circos.trackPlotRegion(
      track.index = 1, 
      panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        
        # Add system labels
        circos.text(
          mean(xlim), ylim[1] + 0.1, 
          sector.name, 
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0, 0.5),
          cex = 1,
          font = 4  # Italic font for consistency
        )
      }, 
      bg.border = NA
    )
    
    # Add title
    title(paste(tool_name, "Defence System Co-occurrence Network"), 
          line = -0.5,  
          cex.main = 1, 
          font.main = 2,
          adj = 0)
    
    # Add legend directly
    legend(x = 1.02, y = -0.1, 
           legend = names(system_colors_subset),
           fill = system_colors_subset,
           title = "Defence Systems",
           cex = 0.8,
           xpd = TRUE
    )
    # Close plot device
    dev.off()
    
    # ========== PART 2: STATISTICAL CORRELATION MATRIX ==========
    
    # Initialize matrices for Fisher's exact test
    odds_ratio_matrix <- matrix(NA, nrow = length(top_systems), ncol = length(top_systems))
    pval_matrix <- matrix(NA, nrow = length(top_systems), ncol = length(top_systems))
    rownames(odds_ratio_matrix) <- colnames(odds_ratio_matrix) <- top_systems
    rownames(pval_matrix) <- colnames(pval_matrix) <- top_systems
    
    # Get total number of genomes
    total_genomes <- nrow(presence_matrix)
    
    # Perform Fisher's exact test for each pair
    for (i in 1:length(top_systems)) {
      for (j in 1:length(top_systems)) {
        # For diagonal (same system), set special values
        if (i == j) {
          odds_ratio_matrix[i, j] <- 100  # Very strong self-association
          pval_matrix[i, j] <- 0.0001     # Highly significant
        } else {
          sys1 <- top_systems[i]
          sys2 <- top_systems[j]
          
          # Create 2x2 contingency table
          # [sys1 present, sys2 present] | [sys1 present, sys2 absent]
          # [sys1 absent, sys2 present]  | [sys1 absent, sys2 absent]
          n11 <- sum(presence_matrix[, sys1] == 1 & presence_matrix[, sys2] == 1)  # Both present
          n10 <- sum(presence_matrix[, sys1] == 1 & presence_matrix[, sys2] == 0)  # sys1 present, sys2 absent
          n01 <- sum(presence_matrix[, sys1] == 0 & presence_matrix[, sys2] == 1)  # sys1 absent, sys2 present
          n00 <- sum(presence_matrix[, sys1] == 0 & presence_matrix[, sys2] == 0)  # Both absent
          
          # Create contingency table
          contingency_table <- matrix(c(n11, n01, n10, n00), nrow = 2, byrow = TRUE)
          
          # Perform Fisher's exact test
          test_result <- fisher.test(contingency_table)
          
          # Store results
          odds_ratio_matrix[i, j] <- test_result$estimate
          pval_matrix[i, j] <- test_result$p.value
        }
      }
    }
    
    # Apply FDR correction
    pval_adjusted_matrix <- matrix(p.adjust(pval_matrix, method = "BH"), 
                                   nrow = length(top_systems), 
                                   ncol = length(top_systems))
    rownames(pval_adjusted_matrix) <- colnames(pval_adjusted_matrix) <- top_systems
    
    # Calculate log2 odds ratio (with limits for infinite values)
    log2_odds_ratio_matrix <- log2(odds_ratio_matrix)
    log2_odds_ratio_matrix[is.infinite(log2_odds_ratio_matrix) & log2_odds_ratio_matrix > 0] <- 16
    log2_odds_ratio_matrix[is.infinite(log2_odds_ratio_matrix) & log2_odds_ratio_matrix < 0] <- -16
    log2_odds_ratio_matrix[is.na(log2_odds_ratio_matrix)] <- 0
    
    # Convert matrices to data frames for ggplot
    fisher_df <- data.frame(
      System1 = rep(top_systems, each = length(top_systems)),
      System2 = rep(top_systems, times = length(top_systems)),
      OddsRatio = as.vector(odds_ratio_matrix),
      Log2OddsRatio = as.vector(log2_odds_ratio_matrix),
      Pvalue = as.vector(pval_matrix),
      Padj = as.vector(pval_adjusted_matrix)
    )
    
    # Set significance threshold
    sig_level <- 0.05
    
    # Create correlation matrix plot
    matrix_plot <- ggplot(fisher_df, aes(x = factor(System1, levels = top_systems), 
                                         y = factor(System2, levels = rev(top_systems)))) +
      # Add background for significant correlations
      geom_tile(
        aes(fill = Padj < sig_level),
        alpha = 0.2,
        color = "grey90"
      ) +
      scale_fill_manual(values = c("white", "yellow"), guide = "none") +
      # Add circles for odds ratio values
      geom_point(
        aes(size = -log10(Padj), color = Log2OddsRatio),
        shape = 16  # Solid circle
      ) +
      # Scale settings
      scale_color_distiller(
        palette = "RdBu",     # ColorBrewer's CVD-safe diverging
        direction = 1,
        limits = c(-16, 16),
        name = expression(log[2](OR))  # Proper mathematical notation
      ) +
      scale_size_continuous(
        range = c(1, 10),
        name = expression(-log[10](p)),  # Proper mathematical notation
        breaks = c(1, 2.5, 5, 10)
      ) +
      # Appearance settings
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title = element_blank(),
        panel.grid = element_line(color = "gray95"),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        legend.key = element_rect(fill = "white", color = NA)
      ) +
      labs(
        title = paste("B.", tool_name, "Defence System Co-occurrence Matrix"),
        subtitle = "Yellow background indicates statistical significance (p < 0.05, FDR-corrected)"
      )
    
    # ========== PART 3: SAVE INDIVIDUAL AND COMBINED PLOTS ==========
    
    # Load circos plot as a grob
    circos_grob <- rasterGrob(readPNG(circos_file), interpolate = TRUE)
    
    # Save individual circos plot
    individual_circos_file <- file.path(output_dir, paste0(output_prefix, "_circos.png"))
    png(individual_circos_file, width = 10, height = 8, units = "in", res = 300)
    grid.draw(circos_grob)
    dev.off()
    
    # Save individual matrix plot
    ggsave(file.path(output_dir, paste0(output_prefix, "_matrix.png")),
           matrix_plot, width = 10, height = 8, dpi = 300)
    
    # Save PDF versions
    ggsave(file.path(output_dir, paste0(output_prefix, "_matrix.pdf")),
           matrix_plot, width = 10, height = 8)
    
    # Combine plots (circos and matrix for this tool)
    combined_plot <- wrap_elements(full = circos_grob) / matrix_plot +
      plot_layout(heights = c(1.2, 1)) +
      plot_annotation(
        title = paste(tool_name, "Defence System Co-occurrence Analysis"),
        theme = theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10)
        )
      )
    
    # Save combined plot (just this tool)
    ggsave(file.path(output_dir, paste0(output_prefix, "_combined.png")),
           combined_plot, width = 10, height = 14, dpi = 300)
    ggsave(file.path(output_dir, paste0(output_prefix, "_combined.pdf")),
           combined_plot, width = 10, height = 14)
    
    # Return important objects
    return(list(
      presence_matrix = presence_matrix,
      cooccur_counts = cooccur_counts,
      odds_ratio_matrix = odds_ratio_matrix,
      pval_adjusted_matrix = pval_adjusted_matrix,
      circos_grob = circos_grob,
      matrix_plot = matrix_plot,
      combined_plot = combined_plot
    ))
  }
  
  # Run co-occurrence analysis for both tools
   defense_cooccur <- analyze_cooccurrence(
    data = defense_df_clean,
    system_col = "type",
    tool_name = "DefenseFinder",
    output_prefix = "defensefinder_cooccur"
  )
  
  padloc_cooccur <- analyze_cooccurrence(
    data = padloc_df_simplified,
    system_col = "Simple_System_Type",
    tool_name = "PADLOC",
    output_prefix = "padloc_cooccur"
  )
  
  # Clean up temporary files
  file.remove(list.files(temp_dir, full.names = TRUE))