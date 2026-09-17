# =============================================================================
# Local copy of the continuous-interpolation helpers used in
# R/cartesian_downsampling/07_Visualise_Continuous_Grid.R
#
# Duplicated (not sourced) so that the composer in this directory stays
# self-contained and does not risk executing 07_'s top-level pipeline code by
# accident. Keep in sync with 07_ if the algorithm changes.
# =============================================================================

# Fill NaN cells inside an interpolated surface from nearest non-NaN neighbour.
fill_interp_holes <- function(z) {
  if (is.null(z)) return(z)
  if (!any(is.na(z))) return(z)
  na_mask <- is.na(z)
  ok_idx  <- which(!na_mask, arr.ind = TRUE)
  if (nrow(ok_idx) == 0) return(z)
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

# Bilinear interpolation on a complete rectangular grid.
bilinear_interp_grid <- function(x_grid, y_grid, z_mat, xo, yo) {
  stopifnot(length(x_grid) == nrow(z_mat),
            length(y_grid) == ncol(z_mat))
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

# Frontier concave hull around measured points.
frontier_concave_hull <- function(x_measured, y_measured, concavity = 0.5) {
  n <- length(x_measured)
  if (n < 3) return(list(hull_x = x_measured, hull_y = y_measured))

  unique_x <- sort(unique(x_measured))
  frontier_points <- list()
  for (x_val in unique_x) {
    at_x <- which(abs(x_measured - x_val) < 1e-6)
    if (length(at_x) > 0) {
      max_y_idx <- at_x[which.max(y_measured[at_x])]
      frontier_points[[length(frontier_points) + 1]] <- max_y_idx
    }
  }

  if (length(frontier_points) < 3) {
    hull_idx <- grDevices::chull(x_measured, y_measured)
    return(list(
      hull_x = c(x_measured[hull_idx], x_measured[hull_idx][1]),
      hull_y = c(y_measured[hull_idx], y_measured[hull_idx][1])
    ))
  }

  frontier_idx <- unlist(frontier_points)
  frontier_x <- x_measured[frontier_idx]
  frontier_y <- y_measured[frontier_idx]
  ord <- order(frontier_x)
  frontier_x <- frontier_x[ord]
  frontier_y <- frontier_y[ord]

  simplify_concave <- function(x_pts, y_pts, concavity) {
    n <- length(x_pts)
    if (n <= 3) return(list(x = x_pts, y = y_pts))
    x_range <- max(x_pts) - min(x_pts)
    y_range <- max(y_pts) - min(y_pts)
    max_edge <- concavity * sqrt(x_range^2 + y_range^2)
    keep <- rep(TRUE, n)
    changed <- TRUE
    iterations <- 0
    while (changed && iterations < 100) {
      changed <- FALSE
      iterations <- iterations + 1
      for (i in 2:(n - 1)) {
        if (!keep[i]) next
        left  <- max(which(keep[1:(i - 1)]))
        right <- i + min(which(keep[(i + 1):n]))
        if (right > n) next
        direct_dist <- sqrt((x_pts[right] - x_pts[left])^2 +
                            (y_pts[right] - y_pts[left])^2)
        path_dist <- sqrt((x_pts[i] - x_pts[left])^2 +
                          (y_pts[i] - y_pts[left])^2) +
                     sqrt((x_pts[right] - x_pts[i])^2 +
                          (y_pts[right] - y_pts[i])^2)
        if (direct_dist < concavity * path_dist && direct_dist < max_edge) {
          keep[i] <- FALSE
          changed <- TRUE
        }
      }
    }
    list(x = x_pts[keep], y = y_pts[keep])
  }

  simplified <- simplify_concave(frontier_x, frontier_y, concavity)
  min_y <- min(y_measured)
  max_x <- max(x_measured)
  min_x <- min(x_measured)
  hull_x <- c(simplified$x, max_x, max_x, min_x, simplified$x[1])
  hull_y <- c(simplified$y, min_y,
              min_y - 0.1 * (max(simplified$y) - min_y),
              min_y - 0.1 * (max(simplified$y) - min_y),
              simplified$y[1])
  list(hull_x = hull_x, hull_y = hull_y)
}

# Test whether interpolation grid points fall inside the concave frontier hull.
points_in_frontier_hull <- function(x_query, y_query, x_measured, y_measured) {
  hull <- frontier_concave_hull(x_measured, y_measured)
  if (requireNamespace("sp", quietly = TRUE)) {
    return(sp::point.in.polygon(x_query, y_query,
                                hull$hull_x, hull$hull_y) > 0)
  }
  n_hull <- length(hull$hull_x)
  result <- rep(FALSE, length(x_query))
  for (i in seq_along(x_query)) {
    x <- x_query[i]; y <- y_query[i]
    inside <- FALSE
    j <- n_hull
    for (k in seq_len(n_hull)) {
      xi <- hull$hull_x[k]; yi <- hull$hull_y[k]
      xj <- hull$hull_x[j]; yj <- hull$hull_y[j]
      if (((yi > y) != (yj > y)) &&
          (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
        inside <- !inside
      }
      j <- k
    }
    result[i] <- inside
  }
  result
}

# Interpolate a cross-score panel (n_donors, y_value, value) onto a dense grid.
# Returns data.frame with columns Donors, Y_Value, .value; values outside the
# frontier concave hull are NA (rendered grey).
interpolate_crossscore_panel <- function(sub, val_col = "value",
                                         n_interp = 200) {
  # Force to plain data.frame so downstream base R subsetting is unambiguous.
  sub <- as.data.frame(sub)
  d <- data.frame(donors_num = sub[["n_donors"]],
                  y_num      = sub[["y_value"]],
                  val_num    = sub[[val_col]])
  d <- d[stats::complete.cases(d[, c("donors_num", "y_num")]), , drop = FALSE]
  measured <- d[!is.na(d$val_num), , drop = FALSE]
  if (nrow(measured) < 4) return(NULL)

  full_grid <- expand.grid(
    donors_num = sort(unique(d$donors_num)),
    y_num      = sort(unique(d$y_num)),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  full_grid <- merge(full_grid, measured, by = c("donors_num", "y_num"),
                     all.x = TRUE, sort = FALSE)

  donors_range <- range(full_grid$donors_num)
  y_range      <- range(full_grid$y_num)
  donors_grid  <- seq(donors_range[1], donors_range[2], length.out = n_interp)
  y_grid       <- seq(y_range[1], y_range[2], length.out = n_interp)

  x_unique <- sort(unique(full_grid$donors_num))
  y_unique <- sort(unique(full_grid$y_num))
  z_mat <- matrix(NA_real_, nrow = length(x_unique), ncol = length(y_unique))
  ix <- match(full_grid$donors_num, x_unique)
  iy <- match(full_grid$y_num,      y_unique)
  for (k in seq_len(nrow(full_grid))) z_mat[ix[k], iy[k]] <- full_grid$val_num[k]

  interp <- bilinear_interp_grid(x_unique, y_unique, z_mat,
                                 donors_grid, y_grid)
  out <- expand.grid(Donors = interp$x, Y_Value = interp$y,
                     KEEP.OUT.ATTRS = FALSE)
  out$.value <- as.vector(interp$z)
  inside <- points_in_frontier_hull(out$Donors, out$Y_Value,
                                    measured$donors_num, measured$y_num)
  out$.value[!inside] <- NA
  out
}

# Scientific power labels for reads-mode Y axis (matches 07_).
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
