# ==============================================================================
# SCRIPT: prettify_labels.R
# PURPOSE: Shared utility for prettifying labels across all figure-generating
#          scripts. Converts internal codes to human-readable labels.
# USAGE: source("R/shared/prettify_labels.R")
# ==============================================================================

#' Prettify a single label or vector of labels
#'
#' Converts internal codes (cell types, datasets, modes, etc.) to human-readable
#' labels suitable for publication figures.
#'
#' @param x Character vector of labels to prettify
#' @return Character vector of prettified labels
#' @examples
#' prettify_label("Astro")        # "Astrocytes"
#' prettify_label("Micro")        # "Microglia"
#' prettify_label("Ground_Truth") # "Ground Truth"
#' prettify_label(c("Astro", "Oligo")) # c("Astrocytes", "Oligodendrocytes")
prettify_label <- function(x) {
  # Define mapping table
  mapping <- c(
    # Cell types
    "Astro" = "Astrocytes",
    "Micro" = "Microglia",
    "Oligo" = "Oligodendrocytes",
    
    # Depth-equalized variants
    "Astro_Shallow" = "Astro (Depth-Eq.)",
    "Astro_DepthEq" = "Astro (Depth-Eq.)",
    "Astro_Shallow_Fewer" = "Astro (Depth+Cell-Eq.)",
    
    # Datasets
    "Roussos" = "PsychAD",
    "Roussos_Full" = "PsychAD (Full)",
    "Roussos_Half_A" = "PsychAD Half A",
    "Roussos_Half_B" = "PsychAD Half B",
    "Roussos_DepthEq" = "PsychAD (Depth-Eq.)",
    "Allen" = "SEA-AD",
    "DeJager" = "ROSMAP",
    "Johansen" = "Johansen",
    "Rexach" = "Rexach",
    
    # Run labels (handle pattern-based)
    "Run_1" = "Split 1",
    "Run_2" = "Split 2",
    "Run_3" = "Split 3",
    "Run_4" = "Split 4",
    "Run_5" = "Split 5",
    "Run_6" = "Split 6",
    
    # Modes
    "cells" = "Cells",
    "reads" = "Counts",
    
    # Categories
    "Ground_Truth" = "Ground Truth",
    "Validation" = "Validation",
    "SexChr" = "Sex Chr.",
    "Autosomal" = "Autosomal",
    
    # LFC bins
    "Small" = "Small (1-2)",
    "Medium" = "Medium (2-3)",
    "Large" = "Large (3+)",
    
    # Common metadata
    "sex" = "Sex",
    "age" = "Age",
    "source" = "Source",
    "Source" = "Brain Bank",
    "disease" = "Disease",
    "ethnicity" = "Ethnicity",
    "subtype" = "Subtype",
    "cell_count" = "Cell Count",
    "cell_type" = "Cell Type",
    "n_donors" = "Number of Donors",
    "n_cells" = "Cells per Donor",
    "n_reads" = "UMIs per Cell",
    "Age_Years" = "Age (Years)",
    "self_reported_ethnicity" = "Self-Reported Ethnicity",
    "Log2_Cells" = "Log2(Cells per Donor)",
    
    # Model predictors (log-transformed)
    "log_donors" = "Number of Donors",
    "log_n_cells" = "Cells per Donor",
    "log_reads_per_cell" = "UMIs per Cell",
    "log_abs_gt_lfc" = "|LFC| (Effect Size)",
    "log_gt_basemean" = "Gene Counts",
    "log_gt_dispersion" = "Dispersion",
    
    # Alternative predictor names
    "abs_gt_lfc" = "|LFC| (Effect Size)",
    "gt_basemean" = "Gene Counts",
    "gt_dispersion" = "Dispersion",
    "reads_per_cell" = "UMIs per Cell",
    "abs_lfc" = "|LFC| (Effect Size)",
    "basemean" = "Gene Counts",
    "dispersion" = "Dispersion"
  )
  
  # Apply mapping, falling back to generic prettification
  result <- vapply(x, function(label) {
    if (is.na(label)) return(NA_character_)
    
    # Check direct mapping first
    if (label %in% names(mapping)) {
      return(mapping[[label]])
    }
    
    # Handle pattern-based transformations
    
    # Roussos_Half_A_R1 -> Roussos Half A (Run 1)
    if (grepl("^Roussos_Half_[AB]_R[0-9]+$", label)) {
      half <- sub(".*Half_([AB]).*", "\\1", label)
      run <- sub(".*_R([0-9]+)$", "\\1", label)
      return(sprintf("PsychAD Half %s (Split %s)", half, run))
    }
    
    # Run_N patterns not in direct map
    if (grepl("^Run_[0-9]+$", label)) {
      run <- sub("Run_", "", label)
      return(sprintf("Split %s", run))
    }
    
    # Generic: replace underscores with spaces
    prettified <- gsub("_", " ", label)
    
    # Title case for single words that weren't mapped
    if (!grepl(" ", prettified)) {
      prettified <- paste0(toupper(substr(prettified, 1, 1)), 
                           substr(prettified, 2, nchar(prettified)))
    }
    
    return(prettified)
  }, character(1), USE.NAMES = FALSE)
  
  return(result)
}

#' Create a labeller function for ggplot2 facets
#'
#' Returns a labeller function compatible with facet_wrap/facet_grid
#'
#' @return A labeller function
#' @examples
#' ggplot(df, aes(x, y)) + 
#'   geom_point() + 
#'   facet_wrap(~cell_type, labeller = prettify_labeller())
prettify_labeller <- function() {
  function(labels) {
    lapply(labels, prettify_label)
  }
}

#' Prettify factor levels in-place
#'
#' Converts factor levels to prettified versions while preserving order
#'
#' @param x Factor to prettify
#' @return Factor with prettified levels
#' @examples
#' df$cell_type <- prettify_factor(df$cell_type)
prettify_factor <- function(x) {
  if (!is.factor(x)) {
    x <- factor(x)
  }
  levels(x) <- prettify_label(levels(x))
  return(x)
}

#' Prettify metadata data values
#'
#' Remaps raw metadata values from datasets to publication-ready labels.
#' Handles disease ontology terms, brain bank codes, ethnicity labels, etc.
#'
#' @param x Character vector of data values to prettify
#' @param column Optional column name hint for context-specific mappings
#' @return Character vector of prettified values
prettify_metadata_values <- function(x, column = NULL) {
  # Disease ontology terms -> readable labels
  disease_map <- c(
    "normal" = "Control",
    "central nervous system disorder" = "CNS Disorder",
    "Alzheimer disease" = "Alzheimer's Disease",
    "Alzheimer's disease" = "Alzheimer's Disease",
    "dementia" = "Dementia"
  )

  # Brain bank source codes -> full names
  source_map <- c(
    "M" = "MSSM",
    "H" = "HBCC",
    "R" = "RADC"
  )

  # Ethnicity labels -> cleaner versions
  ethnicity_map <- c(
    "European" = "European",
    "African" = "African",
    "East Asian" = "East Asian",
    "South Asian" = "South Asian",
    "Asian" = "Asian",
    "Hispanic or Latin American" = "Hispanic/Latino",
    "unknown" = "Unknown",
    "multiethnic" = "Multi-ethnic"
  )

  # Apply context-specific mapping if column is known
  if (!is.null(column)) {
    map <- switch(column,
      "disease" = disease_map,
      "self_reported_ethnicity" = ethnicity_map,
      "Source" = source_map,
      NULL
    )
    if (!is.null(map)) {
      return(ifelse(x %in% names(map), map[x], x))
    }
  }

  # Fall back to trying all maps
  all_maps <- c(disease_map, ethnicity_map, source_map)
  ifelse(x %in% names(all_maps), all_maps[x], x)
}

#' Get standard cell type colour palette
#'
#' Returns a named vector of colours for the three main cell types
#'
#' @param prettified Logical; if TRUE, returns names as prettified labels
#' @return Named character vector of hex colours
#' @examples
#' scale_colour_manual(values = celltype_colours())
#' scale_fill_manual(values = celltype_colours(prettified = TRUE))
celltype_colours <- function(prettified = FALSE) {
  colours <- c(
    "Astro" = "#E69F00",   # Orange
    "Micro" = "#56B4E9",   # Sky blue
    "Oligo" = "#009E73"    # Teal green
  )
  
  if (prettified) {
    names(colours) <- prettify_label(names(colours))
  }
  
  return(colours)
}

#' Get standard cell type colour palette including depth-equalized
#'
#' Returns a named vector of colours including Astro_Shallow
#' 
#' @param prettified Logical; if TRUE, returns names as prettified labels
#' @return Named character vector of hex colours
celltype_colours_extended <- function(prettified = FALSE) {
  colours <- c(
    "Astro" = "#E69F00",         # Orange
    "Micro" = "#56B4E9",         # Sky blue
    "Oligo" = "#009E73",         # Teal green
    "Astro_Shallow" = "#D55E00"  # Red-orange (darker orange)
  )
  
  if (prettified) {
    names(colours) <- prettify_label(names(colours))
  }
  
  return(colours)
}

#' Standard base theme for publication figures
#'
#' Returns a ggplot2 theme with base_size = 16 suitable for paper figures
#' rendered at 32% width
#'
#' @param base_size Base font size (default 16)
#' @param ... Additional arguments passed to theme_bw()
#' @return A ggplot2 theme object
#' @examples
#' ggplot(df, aes(x, y)) + geom_point() + theme_publication()
theme_publication <- function(base_size = 16, ...) {
  theme_bw(base_size = base_size, ...) +
    theme(
      # Clean up panel
      panel.grid.minor = element_blank(),

      # Ensure legend is readable
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),

      # Axis labels
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 2),

      # Facet strips — clean white with subtle bottom border
      strip.text = element_text(size = base_size - 1, face = "bold"),
      strip.background = element_rect(fill = "white", colour = NA),

      # Title and subtitle
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0, colour = "grey40"),

      # Margins
      plot.margin = margin(10, 15, 10, 10)
    )
}

# Message on source
message("[prettify_labels.R] Loaded prettify_label(), prettify_labeller(), ",
        "prettify_factor(), prettify_metadata_values(), celltype_colours(), ",
        "theme_publication()")
