local({
  # Try relative path first (normal host / local R use), then fall back to the
  # absolute path baked into the Docker image by `COPY renv/activate.R` at build
  # time. This lets the same .Rprofile work in both contexts without rebuilding.
  for (p in c("renv/activate.R", "/home/rstudio/renv/activate.R")) {
    if (file.exists(p)) { source(p); break }
  }
})
