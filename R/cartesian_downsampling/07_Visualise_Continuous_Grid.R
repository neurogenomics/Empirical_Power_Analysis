# ==============================================================================
# SCRIPT: Visualize Continuous Grid (Interpolated Heatmaps)
# PURPOSE: Generate continuous heatmaps with bilinear interpolation between
#          measured grid points to reflect true numeric spacing
# INPUT: Aggregated summary data from Step 3 (Aggregated_Summary_Data.qs)
# OUTPUT: Continuous interpolated heatmaps in Continuous_Viz subdirectory
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  library(akima)  # For bilinear interpolation
  
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Package 'ggnewscale' is required. Install: install.packages('ggnewscale')")
  }
  library(ggnewscale)
  
  if (!requireNamespace("shadowtext", quietly = TRUE)) {
    stop("Package 'shadowtext' is required. Install: install.packages('shadowtext')")
  }
  library(shadowtext)
})

# ==============================================================================
# Helper: fill NaN cells inside an interpolated surface from nearest non-NaN
# neighbour. akima::interp's Delaunay triangulation occasionally leaves NaN
# triangles inside the data hull when grid points are irregularly spaced
# (e.g., the log-cells axis 5,10,20,...,3000). This produces empty triangles
# in the rendered heatmap. We detect NaN cells whose row/column has measured
# data on both sides and fill them from the nearest non-NaN interpolated cell.
fill_interp_holes <- function(z) {
  if (is.null(z)) return(z)
  if (!any(is.na(z))) return(z)
  nx <- nrow(z); ny <- ncol(z)
  na_mask <- is.na(z)
  ok_idx  <- which(!na_mask, arr.ind = TRUE)
  if (nrow(ok_idx) == 0) return(z)

  # Bounding box of non-NA cells defines the data envelope. Any NaN inside this
  # rectangle is an akima Delaunay artefact (a triangular hole inside the hull
  # left when grid points are irregularly spaced) and gets filled from its
  # nearest non-NaN neighbour in (row, col) space. The downstream frontier
  # concave-hull mask still carves out anything genuinely outside the measured
  # data region, so over-filling here is safe — extras get masked out later.
  ix_range <- range(ok_idx[, 1])
  iy_range <- range(ok_idx[, 2])
  z_filled <- z
  for (i in ix_range[1]:ix_range[2]) {
    for (j in iy_range[1]:iy_range[2]) {
      if (!na_mask[i, j]) next
      d <- (ok_idx[, 1] - i) ^ 2 + (ok_idx[, 2] - j) ^ 2
      k <- which.min(d)
      z_filled[i, j] <- z[ok_idx[k, 1], ok_idx[k, 2]]
    }
  }
  z_filled
}

# ==============================================================================
# Helper: bilinear interpolation on a complete rectangular grid.
# Replaces akima::interp (Delaunay) which produces triangular holes when grid
# points are irregularly spaced. Inputs:
#   x_grid, y_grid : sorted unique grid coords (length nx, ny)
#   z_mat          : nx*ny matrix of values; MUST be fully populated
#   xo, yo         : output sample coords (any length)
# Returns a list with x=xo, y=yo, z = length(xo)-by-length(yo) matrix.
# Because z_mat is fully populated, output z has no NA cells - by construction
# bilinear interpolation between any 4 valid corners produces a valid result.
bilinear_interp_grid <- function(x_grid, y_grid, z_mat, xo, yo) {
  stopifnot(length(x_grid) == nrow(z_mat),
            length(y_grid) == ncol(z_mat))
  # NA propagation: if any of the 4 surrounding discrete grid corners is NA,
  # the interpolated cell is also NA (renders as grey via na.value in the
  # downstream geom_raster/scale_fill). This makes the continuous surface
  # show exactly the discrete grid's data envelope: regions where there were
  # no measured values stay grey instead of being NN-extrapolated.
  nx <- length(x_grid); ny <- length(y_grid)
  out <- matrix(NA_real_, nrow = length(xo), ncol = length(yo))
  for (i in seq_along(xo)) {
    x <- xo[i]
    if (x < x_grid[1] || x > x_grid[nx]) next
    ix <- findInterval(x, x_grid, all.inside = TRUE)
    x_lo <- x_grid[ix]; x_hi <- x_grid[ix + 1]
    tx <- if (x_hi == x_lo) 0 else (x - x_lo) / (x_hi - x_lo)
    for (j in seq_along(yo)) {
      y <- yo[j]
      if (y < y_grid[1] || y > y_grid[ny]) next
      iy <- findInterval(y, y_grid, all.inside = TRUE)
      y_lo <- y_grid[iy]; y_hi <- y_grid[iy + 1]
      ty <- if (y_hi == y_lo) 0 else (y - y_lo) / (y_hi - y_lo)
      v00 <- z_mat[ix,     iy    ]
      v10 <- z_mat[ix + 1, iy    ]
      v01 <- z_mat[ix,     iy + 1]
      v11 <- z_mat[ix + 1, iy + 1]
      if (is.na(v00) || is.na(v10) || is.na(v01) || is.na(v11)) next
      out[i, j] <- (1 - tx) * (1 - ty) * v00 +
                   tx       * (1 - ty) * v10 +
                   (1 - tx) * ty       * v01 +
                   tx       * ty       * v11
    }
  }
  list(x = xo, y = yo, z = out)
}

# ==============================================================================
# 1. SETUP
# ==============================================================================

# Source config (uses hardcoded PROJECT_ROOT)
if (!exists("PROJECT_ROOT")) {
  possible_config_paths <- c(
    "R/cartesian_downsampling/config.R",
    "config.R",
    "../config.R"
  )
  
  config_loaded <- FALSE
  for (config_path in possible_config_paths) {
    if (file.exists(config_path)) {
      source(config_path)
      config_loaded <- TRUE
      break
    }
  }
  
  if (!config_loaded) {
    stop("Could not find config.R. Ensure working directory is project root or script directory.")
  }
}

source(file.path(SCRIPT_DIR, "utils.R"))
source(file.path(PROJECT_ROOT, "R", "shared", "prettify_labels.R"))

log_section("STEP 7: Visualize Continuous Grid (Interpolated)")

# Mode-aware configuration
y_col <- get_y_col_name()
y_label <- if (DOWNSAMPLE_MODE == "reads") "Target Counts per Donor" else "Cells per Donor"
log_info("Mode: %s | Y-column: %s", DOWNSAMPLE_MODE, y_col)

# Custom label formatter: use 10^n notation for values >= 1e5 (reads mode only)
scientific_power_labels <- function(x) {
  labels <- sapply(x, function(val) {
    if (is.na(val)) return("NA")
    if (abs(val) >= 1e5) {
      exp <- floor(log10(abs(val)))
      coeff <- val / 10^exp
      if (abs(coeff - round(coeff)) < 0.01 && round(coeff) == 1) {
        paste0("10^", exp)
      } else {
        paste0(round(coeff, 1), " %*% 10^", exp)
      }
    } else {
      as.character(as.integer(val))
    }
  })
  parse(text = labels)
}
y_labels_fn <- if (DOWNSAMPLE_MODE == "reads") scientific_power_labels else waiver()

# Check for akima package
if (!requireNamespace("akima", quietly = TRUE)) {
  log_error("Package 'akima' is required for bilinear interpolation.")
  log_info("Install with: install.packages('akima')")
  stop("Missing required package: akima")
}

# Log mode configuration
LOG_MODE_ENABLED <- TRUE  # Set to FALSE to skip log versions
log_info("Log-log continuous plots enabled: %s", LOG_MODE_ENABLED)

#' Format LFC bin label for publication facet strips
#'
#' Converts internal bin names like "Large (3.0+)" to "|LFC| >= 3.0 (n = 25 genes)"
#' @param bin_name Character, the LFC_Bin name from config (e.g. "Large (3.0+)")
#' @param n_genes Numeric, the gene count for this bin
#' @return Character, formatted label
format_lfc_label <- function(bin_name, n_genes) {
  # Extract numeric bounds from bin name patterns like "Tiny (0.2-1.0)", "Large (3.0+)"
  nums <- regmatches(bin_name, gregexpr("[0-9]+\\.?[0-9]*", bin_name))[[1]]
  if (length(nums) >= 2) {
    lo <- nums[1]; hi <- nums[2]
    sprintf("|LFC| %s\u2013%s", lo, hi)
  } else if (length(nums) == 1 && grepl("\\+", bin_name)) {
    sprintf("|LFC| \u2265 %s", nums[1])
  } else {
    bin_name
  }
}

# ==============================================================================
# 2. INTERPOLATION FUNCTIONS
# ==============================================================================

#' Create frontier-boundary polygon around measured points
#' 
#' Creates a boundary that follows the "leading diagonal" frontier of the data,
#' from the highest cell count (top-left) down to highest donor count (bottom-right).
#' The boundary is concave (curves inward toward the origin), creating a shape
#' where the shortest path between edge points goes ABOVE the curve.
#' This eliminates isolated extreme corners while preserving the main data region.
#' 
#' @param x_measured Vector of measured x coordinates (donors)
#' @param y_measured Vector of measured y coordinates (cells)
#' @param concavity Controls boundary tightness (0-1, higher = tighter)
#' @return List with hull_x, hull_y (boundary polygon vertices in order)
frontier_concave_hull <- function(x_measured, y_measured, concavity = 0.5) {
  
  n <- length(x_measured)
  if (n < 3) {
    # Not enough points for meaningful boundary
    return(list(hull_x = x_measured, hull_y = y_measured))
  }
  
  # Find the frontier points: the "upper envelope" in the diagonal direction
  # We want points that are extreme in the (donors + cells) direction
  # but forming a concave curve from top-left to bottom-right
  
  # Calculate diagonal coordinate (direction perpendicular to leading diagonal)
  # Points with high (y - x) are "above" the diagonal, low (y - x) are "below"
  diagonal_coord <- y_measured - x_measured
  
  # Find frontier by scanning along x (donors) and keeping the point with max y
  # at each x position, then smoothing
  unique_x <- sort(unique(x_measured))
  frontier_points <- list()
  
  for (x_val in unique_x) {
    # Get all points at this x value
    at_x <- which(abs(x_measured - x_val) < 1e-6)
    if (length(at_x) > 0) {
      # Keep the point with maximum y at this x
      max_y_idx <- at_x[which.max(y_measured[at_x])]
      frontier_points[[length(frontier_points) + 1]] <- max_y_idx
    }
  }
  
  if (length(frontier_points) < 3) {
    # Fallback to convex hull if not enough frontier points
    hull_idx <- grDevices::chull(x_measured, y_measured)
    return(list(
      hull_x = c(x_measured[hull_idx], x_measured[hull_idx][1]),
      hull_y = c(y_measured[hull_idx], y_measured[hull_idx][1])
    ))
  }
  
  frontier_idx <- unlist(frontier_points)
  frontier_x <- x_measured[frontier_idx]
  frontier_y <- y_measured[frontier_idx]
  
  # Order frontier points from top-left to bottom-right
  # Sort by x (donors), and for each x, we already have the max y
  frontier_order <- order(frontier_x)
  frontier_x <- frontier_x[frontier_order]
  frontier_y <- frontier_y[frontier_order]
  
  # Make the boundary concave by removing points that create convex (outward) angles
  # We want the boundary to curve toward the origin (concave from top-left view)
  simplify_concave <- function(x_pts, y_pts, concavity) {
    n <- length(x_pts)
    if (n <= 3) return(list(x = x_pts, y = y_pts))
    
    # Calculate max allowed distance based on point spread
    x_range <- max(x_pts) - min(x_pts)
    y_range <- max(y_pts) - min(y_pts)
    max_edge <- concavity * sqrt(x_range^2 + y_range^2)
    
    # Start with all points
    keep <- rep(TRUE, n)
    
    # Iteratively remove points that create long edges (convex sections)
    changed <- TRUE
    iterations <- 0
    while (changed && iterations < 100) {
      changed <- FALSE
      iterations <- iterations + 1
      
      # Check each middle point
      for (i in 2:(n-1)) {
        if (!keep[i]) next
        
        # Find neighbors
        left <- max(which(keep[1:(i-1)]))
        right <- i + min(which(keep[(i+1):n]))
        if (right > n) next
        
        # Check edge length if we removed this point
        direct_dist <- sqrt((x_pts[right] - x_pts[left])^2 + 
                           (y_pts[right] - y_pts[left])^2)
        
        # If direct connection is short enough, remove middle point
        path_dist <- sqrt((x_pts[i] - x_pts[left])^2 + (y_pts[i] - y_pts[left])^2) +
                     sqrt((x_pts[right] - x_pts[i])^2 + (y_pts[right] - y_pts[i])^2)
        
        # Remove if the "shortcut" saves enough distance (concavity threshold)
        if (direct_dist < concavity * path_dist && direct_dist < max_edge) {
          keep[i] <- FALSE
          changed <- TRUE
        }
      }
    }
    
    return(list(x = x_pts[keep], y = y_pts[keep]))
  }
  
  # Simplify to create concave shape
  simplified <- simplify_concave(frontier_x, frontier_y, concavity)
  
  # Add boundary points to close the polygon
  # We want a polygon that includes everything BELOW the frontier
  # So we add: frontier (left to right) → bottom-right corner → bottom-left corner → start
  
  min_y <- min(y_measured)
  max_x <- max(x_measured)
  min_x <- min(x_measured)
  
  # Build polygon: frontier (top) → right edge → bottom → left edge → start
  hull_x <- c(simplified$x, max_x, max_x, min_x, simplified$x[1])
  hull_y <- c(simplified$y, min_y, min_y - 0.1 * (max(simplified$y) - min_y), 
              min_y - 0.1 * (max(simplified$y) - min_y), simplified$y[1])
  
  return(list(hull_x = hull_x, hull_y = hull_y))
}

#' Check if points are inside frontier boundary
#' 
#' Uses the concave frontier hull that hugs the data's leading diagonal.
#' This excludes isolated extreme corners (1500×5 and 600×5) while keeping
#' the main data region for interpolation.
#' 
#' @param x_query Vector of x coordinates to test (interpolated grid)
#' @param y_query Vector of y coordinates to test (interpolated grid)
#' @param x_measured Vector of measured x coordinates (actual data points)
#' @param y_measured Vector of measured y coordinates (actual data points)
#' @return Logical vector: TRUE if inside frontier boundary, FALSE if outside
points_in_frontier_hull <- function(x_query, y_query, x_measured, y_measured) {
  
  # Get frontier-based concave hull
  hull <- frontier_concave_hull(x_measured, y_measured)
  
  if (requireNamespace("sp", quietly = TRUE)) {
    # Use sp::point.in.polygon for robustness
    result <- sp::point.in.polygon(x_query, y_query, hull$hull_x, hull$hull_y) > 0
  } else {
    # Ray casting algorithm for point-in-polygon test
    n_hull <- length(hull$hull_x)
    result <- rep(FALSE, length(x_query))
    
    for (i in seq_along(x_query)) {
      x <- x_query[i]
      y <- y_query[i]
      
      # Ray casting: count crossings with polygon edges
      inside <- FALSE
      j <- n_hull
      
      for (k in 1:n_hull) {
        xi <- hull$hull_x[k]
        yi <- hull$hull_y[k]
        xj <- hull$hull_x[j]
        yj <- hull$hull_y[j]
        
        # Check if edge crosses horizontal ray from point
        if (((yi > y) != (yj > y)) && 
            (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
          inside <- !inside
        }
        j <- k
      }
      
      result[i] <- inside
    }
  }
  
  return(result)
}

#' Perform LOCAL bilinear interpolation on grid data
#' 
#' Interpolates ONLY within local neighborhoods of measured points.
#' Uses distance-based masking: each interpolated point must be within
#' max_dist of at least one measured point. Outside = NA (hard gray boundary).
#' This prevents the "donut artifact" where bilinear interpolation bridges
#' across empty regions connecting distant corners.
#' 
#' @param data Data frame with Donors, y_col, and metric columns
#' @param val_col Column name to interpolate
#' @param n_interp Number of interpolation points per dimension (default 200)
#' @return Data frame with interpolated dense grid (NA outside local neighborhoods)
interpolate_grid <- function(data, val_col, n_interp = 200) {

  # Step 1: Pull out the discrete grid actually present in `data` (ignoring the
  # "Max" row, which is annotated separately).
  grid_present <- data %>%
    filter(!is.na(!!sym(y_col)),
           !is.na(Donors),
           as.character(!!sym(y_col)) != "Max")

  if (nrow(grid_present) == 0) {
    return(NULL)
  }

  grid_present$donors_num <- as.numeric(as.character(grid_present$Donors))
  grid_present$y_num      <- as.numeric(as.character(grid_present[[y_col]]))
  grid_present$val_num    <- as.numeric(grid_present[[val_col]])
  grid_present <- grid_present %>%
    filter(!is.na(donors_num), !is.na(y_num))

  measured <- grid_present %>% filter(!is.na(val_num))
  if (nrow(measured) == 0) {
    return(NULL)
  }

  # Step 2: Construct the FULL Cartesian grid from the unique donor and cell
  # levels actually observed in the data. This catches cases where a (donors,
  # cells) cell is missing entirely from `data` (e.g. because no iteration
  # produced an output there) — without this expansion, our NA back-fill has
  # nothing to operate on, and akima::interp leaves a triangular hole inside
  # the data hull that renders as a white triangle in the heatmap.
  full_grid <- expand.grid(
    donors_num = sort(unique(grid_present$donors_num)),
    y_num      = sort(unique(grid_present$y_num)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # Left-join the measured values onto the full grid.
  full_grid <- merge(
    full_grid,
    measured[, c("donors_num", "y_num", "val_num")],
    by = c("donors_num", "y_num"),
    all.x = TRUE,
    sort  = FALSE
  )

  # Step 3: Leave NA cells as NA. The bilinear interpolator propagates NA
  # whenever any of the 4 surrounding discrete grid corners is NA, so the
  # continuous surface shows exactly the discrete grid's data envelope. Cells
  # outside the measured region (e.g. infeasible high-donor x high-cell
  # corners, or low-power low-donor x low-cell corners with no predicted
  # genes) stay NA and render grey via the downstream scale_fill na.value.
  n_missing <- sum(is.na(full_grid$val_num))
  if (n_missing > 0) {
    log_info("  %d cell(s) missing for %s — kept as NA so they render grey",
             n_missing, val_col)
  }

  grid_full <- full_grid

  # The downstream code (range, akima, frontier) used the variable name
  # `finite_data` — keep that name to minimise diff.
  finite_data <- grid_full
  
  if (nrow(finite_data) < 4) {
    return(NULL)  # Need at least 4 points for interpolation
  }
  
  # Create interpolation grid
  donors_range <- range(finite_data$donors_num)
  y_range <- range(finite_data$y_num)
  
  # Create regular grid for interpolation
  donors_grid <- seq(donors_range[1], donors_range[2], length.out = n_interp)
  y_grid <- seq(y_range[1], y_range[2], length.out = n_interp)
  
  # Bilinear interpolation on the rectangular grid. Inputs are guaranteed to
  # cover every (donor, cell) cell after the NN back-fill above, so no holes
  # are possible.
  tryCatch({
    x_unique <- sort(unique(finite_data$donors_num))
    y_unique <- sort(unique(finite_data$y_num))
    z_mat <- matrix(NA_real_, nrow = length(x_unique), ncol = length(y_unique))
    ix <- match(finite_data$donors_num, x_unique)
    iy <- match(finite_data$y_num, y_unique)
    for (k in seq_len(nrow(finite_data))) {
      z_mat[ix[k], iy[k]] <- finite_data$val_num[k]
    }
    interp_result <- bilinear_interp_grid(x_unique, y_unique, z_mat,
                                          donors_grid, y_grid)

    # Convert to long format data frame
    interp_df <- expand.grid(
      Donors = interp_result$x,
      Y_Value = interp_result$y
    )
    interp_df[[val_col]] <- as.vector(interp_result$z)
    names(interp_df)[names(interp_df) == "Y_Value"] <- y_col
    
    # CRITICAL: Mask points outside frontier boundary
    # Creates a concave boundary that follows the "leading diagonal" frontier
    # from high-cells/low-donors down to low-cells/high-donors. Isolated extreme
    # points (1500×5, 600×5) are excluded - no interpolation islands around them.
    # The boundary curves inward (concave), creating the user's requested shape.
    inside_frontier <- points_in_frontier_hull(
      x_query = interp_df$Donors,
      y_query = interp_df[[y_col]],
      x_measured = finite_data$donors_num,
      y_measured = finite_data$y_num
    )
    
    # Set values outside frontier boundary to NA (hard gray boundary)
    interp_df[[val_col]][!inside_frontier] <- NA
    
    return(interp_df)
    
  }, error = function(e) {
    log_warning("Interpolation failed: %s", e$message)
    return(NULL)
  })
}

#' Get Max values for separate annotation
#' 
#' @param data Data frame with aggregated summary data
#' @param val_col Metric column name
#' @return Data frame with Max values only
extract_max_values <- function(data, val_col) {
  data %>%
    filter(as.character(!!sym(y_col)) == "Max" | is.na(!!sym(y_col))) %>%
    group_by(LFC_Bin, Gene_Category) %>%
    summarise(
      Max_Value = mean(!!sym(val_col), na.rm = TRUE),
      .groups = "drop"
    )
}

filter_contourable_data <- function(data, val_col, x_col = "Donors",
                                    y_col_name = y_col, facet_col = NULL,
                                    threshold = 0.8) {
  needed <- c(x_col, y_col_name, val_col)
  if (!all(needed %in% names(data))) {
    return(data[0, , drop = FALSE])
  }

  is_contourable <- function(df) {
    x <- as.numeric(df[[x_col]])
    y <- as.numeric(df[[y_col_name]])
    z <- as.numeric(df[[val_col]])
    ok <- is.finite(x) & is.finite(y) & is.finite(z)
    if (sum(ok) < 4) return(FALSE)
    x <- x[ok]
    y <- y[ok]
    z <- z[ok]
    dplyr::n_distinct(x) >= 2 &&
      dplyr::n_distinct(y) >= 2 &&
      dplyr::n_distinct(z) >= 2 &&
      min(z) <= threshold &&
      max(z) >= threshold
  }

  if (!is.null(facet_col) && facet_col %in% names(data)) {
    idx_by_facet <- split(seq_len(nrow(data)), data[[facet_col]], drop = TRUE)
    keep <- unlist(lapply(idx_by_facet, function(idx) {
      df <- data[idx, , drop = FALSE]
      if (is_contourable(df)) idx else integer(0)
    }), use.names = FALSE)
    return(data[keep, , drop = FALSE])
  }

  if (is_contourable(data)) data else data[0, , drop = FALSE]
}

# ==============================================================================
# 3. PLOTTING FUNCTIONS
# ==============================================================================

#' Create continuous interpolated heatmap
#' 
#' @param interp_data Interpolated data frame
#' @param raw_data Original measured data points
#' @param val_col Column name for the metric
#' @param title_str Plot title
#' @param sub_str Plot subtitle
#' @param max_data Max values data frame (optional)
#' @param fix_scale Whether to fix color scale to 0-1
#' @return ggplot object
plot_continuous_heatmap <- function(interp_data, raw_data, val_col,
                                       title_str, sub_str, max_data = NULL,
                                       fix_scale = TRUE,
                                       y_col_name = y_col,
                                       y_axis_label = y_label,
                                       legend_title = NULL) {
  
  # Prepare raw data for overlay
  raw_finite <- raw_data %>%
    filter(!is.na(!!sym(y_col)), 
           as.character(!!sym(y_col)) != "Max",
           !is.na(Donors),
           !is.na(!!sym(val_col))) %>%
    mutate(
      donors_num = as.numeric(as.character(Donors)),
      y_num = as.numeric(as.character(!!sym(y_col)))
    ) %>%
    filter(!is.na(donors_num), !is.na(y_num))
  
  # Build facet labels
  bin_labels_df <- raw_data %>%
    group_by(LFC_Bin) %>%
    summarise(
      Mean_GT = round(mean(Ref_GT_Count, na.rm = TRUE), 1),
      Str_GT = dplyr::first(GT_Count_String),
      .groups = "drop"
    )
  
  # Prepare raw data for merging labels (simplified - removed detailed counts)
  raw_plot <- raw_finite %>%
    left_join(bin_labels_df, by = "LFC_Bin") %>%
    mutate(LFC_Display_Label = mapply(format_lfc_label, LFC_Bin, Mean_GT))

  # Order facets
  unique_bins <- rev(names(LFC_INTERVALS))
  ordered_labels <- unique(raw_plot$LFC_Display_Label[order(match(raw_plot$LFC_Bin, unique_bins))])
  raw_plot$LFC_Display_Label <- factor(raw_plot$LFC_Display_Label, levels = ordered_labels)
  
  # Add labels to interpolated data
  # Check if LFC_Bin exists in both dataframes
  if (!"LFC_Bin" %in% names(interp_data)) {
    log_warning("LFC_Bin not found in interp_data, creating from unique bins")
    interp_data$LFC_Bin <- raw_finite$LFC_Bin[1]  # Fallback - should already be set
  }
  
  # Use left_join to match LFC_Bin values - ensure both are character type
  interp_data$LFC_Bin <- as.character(interp_data$LFC_Bin)
  bin_labels_df$LFC_Bin <- as.character(bin_labels_df$LFC_Bin)
  
  interp_labeled <- interp_data %>%
    left_join(bin_labels_df, by = "LFC_Bin")
  
  # Check if join worked
  if (any(is.na(interp_labeled$Mean_GT))) {
    log_warning("Some rows have NA after join - filling with defaults")
    interp_labeled$Mean_GT[is.na(interp_labeled$Mean_GT)] <- 0
    interp_labeled$Str_GT[is.na(interp_labeled$Str_GT)] <- "N/A"
  }
  
  interp_labeled <- interp_labeled %>%
    mutate(LFC_Display_Label = mapply(format_lfc_label, LFC_Bin, Mean_GT))
  interp_labeled$LFC_Display_Label <- factor(interp_labeled$LFC_Display_Label,
                                               levels = ordered_labels)
  
  # Derive legend title from metric name if not provided
  fill_label <- if (!is.null(legend_title)) legend_title else {
    lbl <- sub("^Mean_", "", val_col)
    lbl <- sub("_Total$", "", lbl)
    lbl <- sub("_Harm$", "", lbl)
    gsub("_", " ", lbl)
  }

  contour_data <- filter_contourable_data(
    interp_labeled, val_col,
    y_col_name = y_col_name,
    facet_col = "LFC_Display_Label"
  )

  # Base plot with interpolated surface
  p <- ggplot() +
    geom_raster(data = interp_labeled,
                aes(x = Donors, y = !!sym(y_col), fill = !!sym(val_col))) +
    facet_wrap(~LFC_Display_Label, ncol = 2) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(labels = y_labels_fn) +
    labs(title = title_str, subtitle = sub_str,
         x = "Number of Donors", y = y_axis_label) +
    theme_publication(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      panel.spacing = unit(0.8, "lines"),
      legend.position = "right",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18)
    )

  if (nrow(contour_data) > 0) {
    p <- p +
      geom_contour(data = contour_data,
                   aes(x = Donors, y = !!sym(y_col), z = !!sym(val_col),
                       colour = "0.8 threshold"),
                   breaks = c(0.8),
                   linewidth = 1.0,
                   na.rm = TRUE) +
      scale_colour_manual(name = NULL, values = c("0.8 threshold" = "#00FFFF"),
                          guide = guide_legend(override.aes = list(linewidth = 1.2)))
  } else {
    log_info("Skipping 0.8 contour for %s: no contourable surface", val_col)
  }
  
  # Color scale — grey90 outside interpolation boundary
  if (fix_scale) {
    p <- p + scale_fill_viridis_c(option = "magma", limits = c(0, 1),
                                   na.value = "grey90", name = fill_label)
  } else {
    p <- p + scale_fill_viridis_c(option = "magma", na.value = "grey90", name = fill_label)
  }
  
  # Grid point overlay removed for publication clarity
  
  # Add Max annotation if provided
  if (!is.null(max_data) && nrow(max_data) > 0) {
    # Get max values from raw_data for positioning
    max_donors <- max(raw_finite$donors_num, na.rm = TRUE) * 1.1
    
    # Create max annotation data
    max_annot <- max_data %>%
      left_join(bin_labels_df, by = "LFC_Bin") %>%
      mutate(LFC_Display_Label = mapply(format_lfc_label, LFC_Bin, Mean_GT))
    max_annot$LFC_Display_Label <- factor(max_annot$LFC_Display_Label,
                                          levels = ordered_labels)
    
    # Add max values as text annotation on right side
    # We'll use a secondary approach: annotate in plot area
    max_annot <- max_annot %>%
      mutate(
        x_pos = max_donors,
        y_pos = max(raw_finite$y_num, na.rm = TRUE) * 0.5
      )
    
    # Add vertical strip for Max on the right side
    # This is done by extending the plot and adding a separate geom_tile
    p <- p + 
      geom_shadowtext(data = max_annot,
                      aes(x = x_pos, y = y_pos, label = sprintf("Max:\n%.2f", Max_Value)),
                      color = "black", bg.colour = "white", bg.r = 0.1, size = 3)
  }
  
  return(p)
}

#' Create XIST continuous plot
#' 
#' @param interp_data Interpolated XIST data
#' @param raw_data Original XIST data
#' @return ggplot object
plot_xist_continuous <- function(interp_data, raw_data, y_col_name = y_col, y_axis_label = y_label) {
  
  raw_finite <- raw_data %>%
    filter(!is.na(!!sym(y_col)), 
           as.character(!!sym(y_col)) != "Max",
           !is.na(Donors)) %>%
    mutate(
      donors_num = as.numeric(as.character(Donors)),
      y_num = as.numeric(as.character(!!sym(y_col)))
    ) %>%
    filter(!is.na(donors_num), !is.na(y_num))
  
  # Prepare facet labels
  raw_labeled <- raw_finite %>%
    group_by(LFC_Bin) %>%
    mutate(LFC_Label = format_lfc_label(dplyr::first(LFC_Bin), n())) %>%
    ungroup()
  
  unique_bins <- rev(names(LFC_INTERVALS))
  ordered_labels <- unique(raw_labeled$LFC_Label[order(match(raw_labeled$LFC_Bin, unique_bins))])
  raw_labeled$LFC_Label <- factor(raw_labeled$LFC_Label, levels = ordered_labels)
  
  # Ensure LFC_Bin types match
  interp_data$LFC_Bin <- as.character(interp_data$LFC_Bin)
  
  # Get LFC_Bin labels and observation counts
  lfc_bins_df <- raw_finite %>% 
    distinct(LFC_Bin) %>%
    mutate(LFC_Bin = as.character(LFC_Bin))
  
  n_obs_df <- raw_finite %>% 
    group_by(LFC_Bin) %>% 
    summarise(n_obs = n(), .groups = "drop") %>%
    mutate(LFC_Bin = as.character(LFC_Bin))
  
  # Join with interpolated data
  interp_labeled <- interp_data %>%
    left_join(lfc_bins_df, by = "LFC_Bin") %>%
    left_join(n_obs_df, by = "LFC_Bin") %>%
    mutate(LFC_Label = mapply(format_lfc_label, LFC_Bin, n_obs))
  interp_labeled$LFC_Label <- factor(interp_labeled$LFC_Label, levels = ordered_labels)
  
  p <- ggplot() +
    geom_raster(data = interp_labeled,
                aes(x = Donors, y = !!sym(y_col), fill = Mean_Detection)) +
    facet_wrap(~LFC_Label, ncol = 2) +
    coord_cartesian(clip = "off") +  # Prevent edge clipping
    scale_y_continuous(labels = y_labels_fn) +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1),
                          na.value = "grey90", name = "Detection\nRate") +
    labs(title = paste0("XIST Detection (", GENE_OF_INTEREST, ")"),
         x = "Number of Donors", y = y_axis_label) +
    theme_publication(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      panel.spacing = unit(0.8, "lines"),
      legend.position = "right",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18)
    )
  
  return(p)
}

# ==============================================================================
# 4. PROCESS EACH CELL TYPE
# ==============================================================================

cell_types <- get_cell_types()
log_info("Processing %d cell types: %s", length(cell_types), paste(cell_types, collapse = ", "))

for (cell_type in cell_types) {

  for (method in get_viz_methods()) {

  log_section(sprintf("Cell Type: %s | Method: %s", cell_type, method))

  # Get paths
  paths <- get_paths(cell_type, method = method)

  # Create continuous viz directory
  cont_viz_dir <- file.path(paths$plots_dir, "..", "Continuous_Viz")
  cont_viz_dir <- normalizePath(cont_viz_dir, mustWork = FALSE)

  # Skip if continuous visualisation already complete (no-overwrite mode)
  viz_sentinel <- file.path(cont_viz_dir, ".step7_continuous_viz_complete")
  if (!FORCE_RECALC && file.exists(viz_sentinel)) {
    log_info("Skipping %s — continuous visualisation already complete (sentinel exists)", cell_type)
    next
  }

  if (!dir.exists(cont_viz_dir)) {
    dir.create(cont_viz_dir, recursive = TRUE)
    log_info("Created continuous viz directory: %s", cont_viz_dir)
  }
  
  # Load checkpoint
  checkpoint_file <- file.path(paths$plots_dir, "Aggregated_Summary_Data.qs")
  
  if (!file.exists(checkpoint_file)) {
    log_warning("Checkpoint not found for %s, run Step 3 first", cell_type)
    next
  }
  
  checkpoint_data <- qs::qread(checkpoint_file)
  df_summary <- checkpoint_data$df_summary
  df_xist_summary <- checkpoint_data$df_xist_summary
  
  log_info("Loaded summary data: %d rows", nrow(df_summary))
  
  # ==============================================================================
  # 5. GENERATE CONTINUOUS HEATMAPS
  # ==============================================================================
  
  metrics_to_plot <- list(
    list(col = "Mean_Precision", desc = "Mean Precision", sub = NULL),
    list(col = "Mean_Recall", desc = "Mean Recall", sub = NULL)
  )
  
  for (cat in CHR_CATEGORIES) {
    
    df_cat <- df_summary %>% filter(Gene_Category == cat)
    
    if (nrow(df_cat) == 0) {
      log_warning("No data for category %s", cat)
      next
    }
    
    log_info("Processing category: %s (%d rows)", cat, nrow(df_cat))
    
    for (m in metrics_to_plot) {
      
      val_col <- m$col
      
      # Skip if column doesn't exist or has no data
      if (!val_col %in% names(df_cat)) {
        log_warning("Column %s not found, skipping", val_col)
        next
      }
      
      if (all(is.na(df_cat[[val_col]]))) {
        log_warning("All values NA for %s in %s, skipping", val_col, cat)
        next
      }
      
      # Interpolate for each LFC bin
      interp_results <- list()
      max_values <- list()
      
      for (bin_name in unique(df_cat$LFC_Bin)) {
        
        df_bin <- df_cat %>% filter(LFC_Bin == bin_name)
        
        # Interpolate finite grid
        interp_df <- interpolate_grid(df_bin, val_col, n_interp = 200)
        
        if (!is.null(interp_df)) {
          interp_df$LFC_Bin <- bin_name
          interp_df$Gene_Category <- cat
          interp_results[[length(interp_results) + 1]] <- interp_df
        }
        
        # Extract Max values
        max_df <- extract_max_values(df_bin, val_col)
        if (nrow(max_df) > 0) {
          max_df$LFC_Bin <- bin_name
          max_values[[length(max_values) + 1]] <- max_df
        }
      }
      
      if (length(interp_results) == 0) {
        log_warning("No interpolation results for %s in %s", val_col, cat)
        next
      }
      
      # Combine interpolation results
      interp_combined <- do.call(rbind, interp_results)
      max_combined <- if (length(max_values) > 0) do.call(rbind, max_values) else NULL
      
      # Create plot
      title_str <- paste0(prettify_label(cell_type), ": ", m$desc)

      p <- plot_continuous_heatmap(
        interp_data = interp_combined,
        raw_data = df_cat,
        val_col = val_col,
        title_str = title_str,
        sub_str = m$sub,
        max_data = max_combined,
        fix_scale = TRUE,
        y_col_name = y_col,
        y_axis_label = y_label
      )
      
      # Save
      fname <- paste0(val_col, "_", cat, "_Continuous.png")
      fpath <- file.path(cont_viz_dir, fname)
      
      ggsave(fpath, p, width = 14, height = 10, dpi = 150)
      log_info("Saved: %s", fname)
    }
  }
  
  # ==============================================================================
  # 5b. GENERATE INDIVIDUAL LFC BIN PLOTS (Separate full-page plots per bin)
  # ==============================================================================
  
  log_info("Generating individual LFC bin plots (separate full-page plots)...")
  
  # Create subdirectory for individual bin plots
  individual_bins_dir <- file.path(cont_viz_dir, "Individual_Bins")
  if (!dir.exists(individual_bins_dir)) {
    dir.create(individual_bins_dir, recursive = TRUE)
    log_info("Created individual bins directory: %s", individual_bins_dir)
  }
  
  for (cat in CHR_CATEGORIES) {
    
    df_cat <- df_summary %>% filter(Gene_Category == cat)
    if (nrow(df_cat) == 0) next
    
    for (m in metrics_to_plot) {
      
      val_col <- m$col
      if (!val_col %in% names(df_cat)) next
      if (all(is.na(df_cat[[val_col]]))) next
      
      # Generate separate plot for each LFC bin
      for (bin_name in unique(df_cat$LFC_Bin)) {
        
        df_bin <- df_cat %>% filter(LFC_Bin == bin_name)
        if (nrow(df_bin) == 0) next
        
        # Interpolate for this bin only
        interp_df <- interpolate_grid(df_bin, val_col, n_interp = 300)  # Higher res for single bin
        
        if (is.null(interp_df)) {
          log_warning("Could not interpolate for %s %s %s", cat, val_col, bin_name)
          next
        }
        
        interp_df$LFC_Bin <- bin_name
        interp_df$Gene_Category <- cat
        
        # Get GT count info for title (simplified - no detailed counts)
        bin_info <- df_bin %>%
          summarise(
            Mean_GT = round(mean(Ref_GT_Count, na.rm = TRUE), 1),
            .groups = "drop"
          )
        
        # Prepare raw data for overlay
        raw_plot <- df_bin %>%
          filter(!is.na(!!sym(y_col)), 
                 as.character(!!sym(y_col)) != "Max",
                 !is.na(Donors),
                 !is.na(!!sym(val_col))) %>%
          mutate(
            donors_num = as.numeric(as.character(Donors)),
            y_num = as.numeric(as.character(!!sym(y_col)))
          ) %>%
          filter(!is.na(donors_num), !is.na(y_num))
        
        # Derive legend title
        fill_label_ind <- sub("^Mean_", "", val_col)
        fill_label_ind <- sub("_Total$", "", fill_label_ind)
        fill_label_ind <- sub("_Harm$", "", fill_label_ind)
        fill_label_ind <- gsub("_", " ", fill_label_ind)

        contour_data <- filter_contourable_data(interp_df, val_col)

        # Build single-bin plot (no faceting)
        p_individual <- ggplot() +
          geom_raster(data = interp_df,
                      aes(x = Donors, y = !!sym(y_col), fill = !!sym(val_col))) +
          coord_cartesian(clip = "off") +
          scale_fill_viridis_c(option = "magma", limits = c(0, 1),
                               na.value = "grey90", name = fill_label_ind) +
          labs(
            title = paste0(prettify_label(cell_type), ": ", m$desc),
            subtitle = format_lfc_label(bin_name, bin_info$Mean_GT),
            x = "Number of Donors",
            y = y_label
          ) +
          theme_publication(base_size = 18) +
          theme(
            panel.grid = element_blank(),
            legend.position = "right",
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18)
          )

        if (nrow(contour_data) > 0) {
          p_individual <- p_individual +
            geom_contour(data = contour_data,
                         aes(x = Donors, y = !!sym(y_col), z = !!sym(val_col),
                             colour = "0.8 threshold"),
                         breaks = c(0.8),
                         linewidth = 1.2,
                         na.rm = TRUE) +
            scale_colour_manual(name = NULL, values = c("0.8 threshold" = "#00FFFF"),
                                guide = guide_legend(override.aes = list(linewidth = 1.2)))
        } else {
          log_info("Skipping 0.8 contour for %s %s: no contourable surface",
                   val_col, bin_name)
        }
        
        # Save with larger dimensions for single-bin detail
        bin_suffix <- gsub("[^a-zA-Z0-9_-]", "", gsub(" ", "_", bin_name))  # Clean filename
        fname <- paste0(val_col, "_", cat, "_", bin_suffix, "_Individual.png")
        fpath <- file.path(individual_bins_dir, fname)
        
        ggsave(fpath, p_individual, width = 12, height = 10, dpi = 200)  # Higher DPI for detail
        log_info("Saved individual plot: %s", fname)
      }
    }
  }
  
  log_info("Completed individual LFC bin plots")
  
  # ==============================================================================
  # 6. XIST CONTINUOUS PLOT
   # ==============================================================================
  
  log_info("Generating XIST continuous plot...")
  
  interp_xist_results <- list()
  
  for (bin_name in unique(df_xist_summary$LFC_Bin)) {
    df_bin <- df_xist_summary %>% filter(LFC_Bin == bin_name)
    
    interp_df <- interpolate_grid(df_bin, "Mean_Detection", n_interp = 200)
    
    if (!is.null(interp_df)) {
      interp_df$LFC_Bin <- bin_name
      interp_xist_results[[length(interp_xist_results) + 1]] <- interp_df
    }
  }
  
  if (length(interp_xist_results) > 0) {
    interp_xist_combined <- do.call(rbind, interp_xist_results)
    
    p_xist <- plot_xist_continuous(interp_xist_combined, df_xist_summary, 
                                  y_col_name = y_col, y_axis_label = y_label)
    
    ggsave(file.path(cont_viz_dir, "XIST_Sanity_Check_Continuous.png"), 
           p_xist, width = 14, height = 10, dpi = 150)
    log_info("Saved: XIST_Sanity_Check_Continuous.png")
  } else {
    log_warning("Could not generate XIST continuous plot (insufficient data)")
  }
  
  # ==============================================================================
  # 7. READ DEPTH CONTINUOUS PLOT (cells mode only)
  # ==============================================================================
  
  if (DOWNSAMPLE_MODE == "cells" && !is.null(checkpoint_data$df_reads_summary)) {
    
    log_info("Generating read depth continuous plot...")
    
    df_reads <- checkpoint_data$df_reads_summary
    
    # Interpolate reads
    interp_reads_results <- list()
    
    for (bin_name in c("All")) {  # Single "bin" for reads
      interp_df <- interpolate_grid(df_reads, "Mean_Total_Reads", n_interp = 200)
      
      if (!is.null(interp_df)) {
        interp_df$LFC_Bin <- bin_name
        interp_reads_results[[length(interp_reads_results) + 1]] <- interp_df
      }
    }
    
    if (length(interp_reads_results) > 0) {
      interp_reads_combined <- do.call(rbind, interp_reads_results)
      
      # Get finite raw data
      raw_finite <- df_reads %>%
        filter(!is.na(!!sym(y_col)), 
               as.character(!!sym(y_col)) != "Max",
               !is.na(Donors)) %>%
        mutate(
          donors_num = as.numeric(as.character(Donors)),
          y_num = as.numeric(as.character(!!sym(y_col)))
        ) %>%
        filter(!is.na(donors_num), !is.na(y_num))
      
      # Determine global scale
      global_min <- min(interp_reads_combined$Mean_Total_Reads, na.rm = TRUE)
      global_max <- max(interp_reads_combined$Mean_Total_Reads, na.rm = TRUE)
      
      p_reads <- ggplot() +
        geom_raster(data = interp_reads_combined,
                    aes(x = Donors, y = !!sym(y_col), fill = Mean_Total_Reads)) +
        scale_y_continuous(labels = y_labels_fn) +
        scale_fill_viridis_c(option = "plasma", trans = "log10",
                             limits = c(global_min, global_max),
                             na.value = "grey90", name = "Total Counts") +
        labs(title = paste0(prettify_label(cell_type), ": Total Counts per Matrix"),
             x = "Number of Donors", y = y_label) +
        theme_publication(base_size = 16) +
        theme(
          panel.grid = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18)
        )
      
      ggsave(file.path(cont_viz_dir, "Total_Reads_Per_Matrix_Continuous.png"),
             p_reads, width = 10, height = 8, dpi = 150)
      log_info("Saved: Total_Reads_Per_Matrix_Continuous.png")
    }
  }
  
  # Write sentinel so subsequent runs skip this cell type
  writeLines(as.character(Sys.time()), viz_sentinel)
  log_info("Completed continuous visualization for %s/%s (sentinel written)", cell_type, method)
  }
}

# ==============================================================================
# 8. GENERATE SEMI-LOG VERSIONS (Optional) - Log Y-axis (cells), Linear X-axis (donors)
# ==============================================================================

if (LOG_MODE_ENABLED) {
  log_section("Generating Semi-Log Continuous Grid Visualizations")
  log_info("Creating semi-log versions in Continuous_Viz_Log/ (Y: log scale, X: linear)")
  
  for (cell_type in cell_types) {
    for (method in get_viz_methods()) {

    log_section(sprintf("Cell Type: %s | Method: %s (Log-Log)", cell_type, method))

    # Get paths
    paths <- get_paths(cell_type, method = method)

    # Create log continuous viz directory
    cont_viz_log_dir <- file.path(paths$plots_dir, "..", "Continuous_Viz_Log")
    cont_viz_log_dir <- normalizePath(cont_viz_log_dir, mustWork = FALSE)

    # Skip if semi-log visualisation already complete (no-overwrite mode)
    viz_log_sentinel <- file.path(cont_viz_log_dir, ".step7_semilog_viz_complete")
    if (!FORCE_RECALC && file.exists(viz_log_sentinel)) {
      log_info("Skipping %s — semi-log visualisation already complete (sentinel exists)", cell_type)
      next
    }

    if (!dir.exists(cont_viz_log_dir)) {
      dir.create(cont_viz_log_dir, recursive = TRUE)
      log_info("Created log viz directory: %s", cont_viz_log_dir)
    }
    
    # Load checkpoint
    checkpoint_file <- file.path(paths$plots_dir, "Aggregated_Summary_Data.qs")
    
    if (!file.exists(checkpoint_file)) {
      next
    }
    
    checkpoint_data <- qs::qread(checkpoint_file)
    df_summary <- checkpoint_data$df_summary
    
    # Generate log versions
    for (cat in CHR_CATEGORIES) {
      
      df_cat <- df_summary %>% filter(Gene_Category == cat)
      
      if (nrow(df_cat) == 0) next
      
      log_info("Processing log version category: %s", cat)
      
      # Metrics for log version
      metrics_to_plot <- list(
        list(col = "Mean_Recall", desc = "Mean Recall"),
        list(col = "Mean_Precision", desc = "Mean Precision")
      )
      
      for (m in metrics_to_plot) {
        
        val_col <- m$col
        
        if (!val_col %in% names(df_cat)) next
        if (all(is.na(df_cat[[val_col]]))) next
        
        # Interpolate for each LFC bin in LOG space
        interp_results <- list()
        
        for (bin_name in unique(df_cat$LFC_Bin)) {

          df_bin <- df_cat %>% filter(LFC_Bin == bin_name)

          # Build the FULL discrete grid (keep cells with NA val_col so we can
          # back-fill them from neighbouring measured cells — same rule as the
          # discrete grid plot uses).
          all_grid <- df_bin %>%
            filter(!is.na(!!sym(y_col)),
                   !is.na(Donors),
                   as.character(!!sym(y_col)) != "Max") %>%
            mutate(
              donors_num = as.numeric(as.character(Donors)),
              y_num = as.numeric(as.character(!!sym(y_col)))
            ) %>%
            filter(!is.na(donors_num), !is.na(y_num),
                   donors_num > 0, y_num > 0)  # Must be positive for log

          measured <- all_grid %>% filter(!is.na(!!sym(val_col)))
          if (nrow(measured) == 0) next

          # Build the FULL Cartesian grid from observed donor and cell levels
          # so any (donors, cells) combo missing from `df_bin` gets a row to
          # back-fill. Without this expansion akima leaves triangular holes
          # inside the data hull when the input has missing rows (not just
          # rows with NA values).
          full_grid <- expand.grid(
            donors_num = sort(unique(all_grid$donors_num)),
            y_num      = sort(unique(all_grid$y_num)),
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
          )
          full_grid <- merge(
            full_grid,
            measured[, c("donors_num", "y_num", val_col)],
            by = c("donors_num", "y_num"),
            all.x = TRUE,
            sort  = FALSE
          )

          # Leave NA cells as NA: bilinear_interp_grid will propagate NA so the
          # continuous surface mirrors the discrete grid's data envelope.
          n_missing <- sum(is.na(full_grid[[val_col]]))
          if (n_missing > 0) {
            log_info("  [SemiLog %s/%s] %d cell(s) missing — kept as NA",
                     val_col, bin_name, n_missing)
          }

          finite_data <- full_grid

          if (nrow(finite_data) < 4) next

          # Transform Y to log space, keep X linear (semi-log interpolation)
          log_finite <- finite_data %>%
            mutate(
              # Keep donors linear (no log transform)
              y_log = log10(y_num)
            )

          # Interpolate in SEMI-LOG space (X linear, Y log)
          donors_range <- range(log_finite$donors_num)  # X: linear
          y_range <- range(log_finite$y_log)            # Y: log

          donors_grid <- seq(donors_range[1], donors_range[2], length.out = 200)
          y_grid <- seq(y_range[1], y_range[2], length.out = 200)

          # Bilinear on rectangular grid (no Delaunay, no internal holes).
          interp_result <- tryCatch({
            x_unique <- sort(unique(log_finite$donors_num))
            y_unique <- sort(unique(log_finite$y_log))
            z_mat <- matrix(NA_real_, nrow = length(x_unique), ncol = length(y_unique))
            ix <- match(log_finite$donors_num, x_unique)
            iy <- match(log_finite$y_log, y_unique)
            for (k in seq_len(nrow(log_finite))) {
              z_mat[ix[k], iy[k]] <- log_finite[[val_col]][k]
            }
            bilinear_interp_grid(x_unique, y_unique, z_mat, donors_grid, y_grid)
          }, error = function(e) NULL)

          if (is.null(interp_result)) next

          # Convert Y back to original scale, keep X linear
          interp_df <- expand.grid(
            Donors = interp_result$x,       # X: linear (no conversion)
            Y_Value = 10^interp_result$y    # Y: log back to linear
          )
          names(interp_df)[names(interp_df) == "Y_Value"] <- y_col
          interp_df[[val_col]] <- as.vector(interp_result$z)

          # Apply frontier masking in semi-log space (X linear, Y log)
          hull <- frontier_concave_hull(log_finite$donors_num, log_finite$y_log)
          inside <- tryCatch({
            if (requireNamespace("sp", quietly = TRUE)) {
              sp::point.in.polygon(interp_df$Donors, log10(interp_df[[y_col]]),
                                  hull$hull_x, hull$hull_y) > 0
            } else {
              # Simple bounding box fallback (Y in original scale, so compare to log bounds)
              interp_df$Donors >= min(finite_data$donors_num) &
              interp_df$Donors <= max(finite_data$donors_num) &
              interp_df[[y_col]] >= min(finite_data$y_num) &
              interp_df[[y_col]] <= max(finite_data$y_num)
            }
          }, error = function(e) rep(TRUE, nrow(interp_df)))

          interp_df[[val_col]][!inside] <- NA
          
          interp_df$LFC_Bin <- bin_name
          interp_df$Gene_Category <- cat
          interp_results[[length(interp_results) + 1]] <- interp_df
        }
        
        if (length(interp_results) == 0) next
        
        interp_combined <- do.call(rbind, interp_results)
        
        # Create log-log plot
        # Prepare raw data for overlay (log scale)
        raw_plot <- df_cat %>%
          filter(!is.na(!!sym(y_col)), 
                 as.character(!!sym(y_col)) != "Max",
                 !is.na(Donors),
                 !is.na(!!sym(val_col))) %>%
          mutate(
            donors_num = as.numeric(as.character(Donors)),
            y_num = as.numeric(as.character(!!sym(y_col)))
          ) %>%
          filter(!is.na(donors_num), !is.na(y_num),
                 donors_num > 0, y_num > 0) %>%
          group_by(LFC_Bin) %>%
          mutate(
            Mean_GT = round(mean(Ref_GT_Count, na.rm = TRUE), 1)
          ) %>%
          ungroup() %>%
          mutate(LFC_Display_Label = mapply(format_lfc_label, LFC_Bin, Mean_GT))

        unique_bins <- rev(names(LFC_INTERVALS))
        ordered_labels <- unique(raw_plot$LFC_Display_Label[order(match(raw_plot$LFC_Bin, unique_bins))])
        raw_plot$LFC_Display_Label <- factor(raw_plot$LFC_Display_Label, levels = ordered_labels)

        # Add labels to interpolated
        interp_labeled <- interp_combined %>%
          left_join(
            raw_plot %>% distinct(LFC_Bin, LFC_Display_Label),
            by = "LFC_Bin"
          )

        # Derive legend title from metric name
        fill_label_log <- sub("^Mean_", "", val_col)
        fill_label_log <- sub("_Total$", "", fill_label_log)
        fill_label_log <- sub("_Harm$", "", fill_label_log)
        fill_label_log <- gsub("_", " ", fill_label_log)

        contour_data <- filter_contourable_data(
          interp_labeled, val_col,
          y_col_name = y_col,
          facet_col = "LFC_Display_Label"
        )

        # Build plot with semi-log scale (Y log, X linear)
        p_log <- ggplot() +
          geom_raster(data = interp_labeled,
                      aes(x = Donors, y = !!sym(y_col), fill = !!sym(val_col))) +
          facet_wrap(~LFC_Display_Label, ncol = 2) +
          coord_cartesian(clip = "off") +
          scale_y_log10(labels = y_labels_fn) +
          scale_fill_viridis_c(option = "magma", limits = c(0, 1),
                               na.value = "grey90", name = fill_label_log) +
          labs(title = paste0(prettify_label(cell_type), ": ", m$desc),
               x = "Number of Donors",
               y = y_label) +
          theme_publication(base_size = 16) +
          theme(
            panel.grid = element_blank(),
            strip.text = element_text(face = "bold", size = 13),
            panel.spacing = unit(0.8, "lines"),
            # Harmonised with 04_Visualise_Power.R so grid + surface composites
            # in R/figure_composition have matching legend widths.
            legend.position = "right",
            legend.key.width = unit(0.7, "cm"),
            legend.key.height = unit(1.3, "cm"),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18)
          )

        if (nrow(contour_data) > 0) {
          p_log <- p_log +
            geom_contour(data = contour_data,
                         aes(x = Donors, y = !!sym(y_col), z = !!sym(val_col),
                             colour = "0.8 threshold"),
                         breaks = c(0.8),
                         linewidth = 1.0,
                         na.rm = TRUE) +
            scale_colour_manual(name = NULL, values = c("0.8 threshold" = "#00FFFF"),
                                guide = guide_legend(override.aes = list(linewidth = 1.2)))
        } else {
          log_info("Skipping semi-log 0.8 contour for %s: no contourable surface",
                   val_col)
        }

        # Save
        fname <- paste0(val_col, "_", cat, "_Continuous_SemiLog.png")
        fpath <- file.path(cont_viz_log_dir, fname)
        
        ggsave(fpath, p_log, width = 14, height = 10, dpi = 150)
        log_info("Saved: %s", fname)
      }
    }
    
    log_info("Completed semi-log faceted visualizations for %s/%s", cell_type, method)
    
    # ==============================================================================
    # 8b. GENERATE INDIVIDUAL LFC BIN PLOTS FOR SEMI-LOG VERSION
    # ==============================================================================
    
    log_info("Generating individual LFC bin plots for semi-log version...")
    
    # Create subdirectory for log-log individual bin plots
    log_individual_bins_dir <- file.path(cont_viz_log_dir, "Individual_Bins")
    if (!dir.exists(log_individual_bins_dir)) {
      dir.create(log_individual_bins_dir, recursive = TRUE)
      log_info("Created log individual bins directory: %s", log_individual_bins_dir)
    }
    
    # Only Total category, only Recall and Precision
    log_metrics <- list(
      list(col = "Mean_Recall", desc = "Mean Recall"),
      list(col = "Mean_Precision", desc = "Mean Precision")
    )
    
    for (cat in c("Total")) {
      
      df_cat <- df_summary %>% filter(Gene_Category == cat)
      if (nrow(df_cat) == 0) next
      
      for (m in log_metrics) {
        
        val_col <- m$col
        if (!val_col %in% names(df_cat)) next
        if (all(is.na(df_cat[[val_col]]))) next
        
        # Generate separate plot for each LFC bin
        for (bin_name in unique(df_cat$LFC_Bin)) {

          df_bin <- df_cat %>% filter(LFC_Bin == bin_name)
          if (nrow(df_bin) == 0) next

          # Build the FULL discrete grid (NA val_col cells included for back-fill).
          all_grid <- df_bin %>%
            filter(!is.na(!!sym(y_col)),
                   as.character(!!sym(y_col)) != "Max",
                   !is.na(Donors)) %>%
            mutate(
              donors_num = as.numeric(as.character(Donors)),
              y_num = as.numeric(as.character(!!sym(y_col)))
            ) %>%
            filter(!is.na(donors_num), !is.na(y_num),
                   donors_num > 0, y_num > 0)

          measured <- all_grid %>% filter(!is.na(!!sym(val_col)))
          if (nrow(measured) == 0) next

          # Build the FULL Cartesian grid so any (donors, cells) combo absent
          # from `df_bin` gets a row to back-fill.
          full_grid <- expand.grid(
            donors_num = sort(unique(all_grid$donors_num)),
            y_num      = sort(unique(all_grid$y_num)),
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
          )
          full_grid <- merge(
            full_grid,
            measured[, c("donors_num", "y_num", val_col)],
            by = c("donors_num", "y_num"),
            all.x = TRUE,
            sort  = FALSE
          )

          # Leave NA cells as NA — bilinear_interp_grid propagates NA, so the
          # continuous surface mirrors the discrete grid's data envelope.
          n_missing <- sum(is.na(full_grid[[val_col]]))
          if (n_missing > 0) {
            log_info("  [SemiLog-Indiv %s/%s] %d cell(s) missing — kept as NA",
                     val_col, bin_name, n_missing)
          }

          finite_data <- full_grid

          if (nrow(finite_data) < 4) next

          # Transform Y to log space, keep X linear (semi-log interpolation)
          log_finite <- finite_data %>%
            mutate(
              # Keep donors linear (no log transform)
              y_log = log10(y_num)
            )

          # Interpolate in SEMI-LOG space (X linear, Y log)
          donors_range <- range(log_finite$donors_num)  # X: linear
          y_range <- range(log_finite$y_log)            # Y: log

          donors_grid <- seq(donors_range[1], donors_range[2], length.out = 300)
          y_grid <- seq(y_range[1], y_range[2], length.out = 300)

          # Bilinear on rectangular grid (no Delaunay, no internal holes).
          interp_result <- tryCatch({
            x_unique <- sort(unique(log_finite$donors_num))
            y_unique <- sort(unique(log_finite$y_log))
            z_mat <- matrix(NA_real_, nrow = length(x_unique), ncol = length(y_unique))
            ix <- match(log_finite$donors_num, x_unique)
            iy <- match(log_finite$y_log, y_unique)
            for (k in seq_len(nrow(log_finite))) {
              z_mat[ix[k], iy[k]] <- log_finite[[val_col]][k]
            }
            bilinear_interp_grid(x_unique, y_unique, z_mat, donors_grid, y_grid)
          }, error = function(e) NULL)

          if (is.null(interp_result)) next

          # Convert Y back to original scale, keep X linear
          interp_df <- expand.grid(
            Donors = interp_result$x,       # X: linear (no conversion)
            Y_Value = 10^interp_result$y    # Y: log back to linear
          )
          names(interp_df)[names(interp_df) == "Y_Value"] <- y_col
          interp_df[[val_col]] <- as.vector(interp_result$z)
          
          # Apply frontier masking in semi-log space (X linear, Y log)
          hull <- frontier_concave_hull(log_finite$donors_num, log_finite$y_log)
          inside <- tryCatch({
            if (requireNamespace("sp", quietly = TRUE)) {
              sp::point.in.polygon(interp_df$Donors, log10(interp_df[[y_col]]), 
                                  hull$hull_x, hull$hull_y) > 0
            } else {
              interp_df$Donors >= min(finite_data$donors_num) &
              interp_df$Donors <= max(finite_data$donors_num) &
              interp_df[[y_col]] >= min(finite_data$y_num) &
              interp_df[[y_col]] <= max(finite_data$y_num)
            }
          }, error = function(e) rep(TRUE, nrow(interp_df)))
          
          interp_df[[val_col]][!inside] <- NA
          
          # Get GT count info (simplified)
          bin_info <- df_bin %>%
            summarise(
              Mean_GT = round(mean(Ref_GT_Count, na.rm = TRUE), 1),
              .groups = "drop"
            )
          
          # Prepare raw data for overlay
          raw_plot <- finite_data
          
          # Derive legend title
          fill_label_slog <- sub("^Mean_", "", val_col)
          fill_label_slog <- sub("_Total$", "", fill_label_slog)
          fill_label_slog <- sub("_Harm$", "", fill_label_slog)
          fill_label_slog <- gsub("_", " ", fill_label_slog)

          contour_data <- filter_contourable_data(interp_df, val_col)

          # Build individual semi-log plot (Y log, X linear)
          p_log_individual <- ggplot() +
            geom_raster(data = interp_df,
                        aes(x = Donors, y = !!sym(y_col), fill = !!sym(val_col))) +
            coord_cartesian(clip = "off") +
            scale_y_log10(labels = y_labels_fn) +
            scale_fill_viridis_c(option = "magma", limits = c(0, 1),
                                 na.value = "grey90", name = fill_label_slog) +
            labs(
              title = paste0(prettify_label(cell_type), ": ", m$desc),
              subtitle = format_lfc_label(bin_name, bin_info$Mean_GT),
              x = "Number of Donors",
              y = y_label
            ) +
            theme_publication(base_size = 18) +
            theme(
              panel.grid = element_blank(),
              legend.position = "right",
              legend.text = element_text(size = 16),
              legend.title = element_text(size = 18)
            )

          if (nrow(contour_data) > 0) {
            p_log_individual <- p_log_individual +
              geom_contour(data = contour_data,
                           aes(x = Donors, y = !!sym(y_col), z = !!sym(val_col),
                               colour = "0.8 threshold"),
                           breaks = c(0.8),
                           linewidth = 1.2,
                           na.rm = TRUE) +
              scale_colour_manual(name = NULL, values = c("0.8 threshold" = "#00FFFF"),
                                  guide = guide_legend(override.aes = list(linewidth = 1.2)))
          } else {
            log_info("Skipping semi-log 0.8 contour for %s %s: no contourable surface",
                     val_col, bin_name)
          }
          
          # Save
          bin_suffix <- gsub("[^a-zA-Z0-9_-]", "", gsub(" ", "_", bin_name))
          fname <- paste0(val_col, "_", cat, "_", bin_suffix, "_Individual_SemiLog.png")
          fpath <- file.path(log_individual_bins_dir, fname)
          
          ggsave(fpath, p_log_individual, width = 12, height = 10, dpi = 200)
          log_info("Saved semi-log individual plot: %s", fname)
        }
      }
    }
    
    log_info("Completed semi-log individual visualizations for %s/%s", cell_type, method)

    # Write sentinel so subsequent runs skip this cell type
    writeLines(as.character(Sys.time()), viz_log_sentinel)
    log_info("Semi-log visualisation sentinel written for %s/%s", cell_type, method)
    }
  }
}

log_section("STEP 7 COMPLETE")
log_info("All continuous visualizations generated in Continuous_Viz/ subdirectories")
if (LOG_MODE_ENABLED) {
  log_info("Semi-log versions (log Y-axis only) available in Continuous_Viz_Log/ subdirectories")
}
