local({
  # Inside Docker (EPA_DOCKER=1), packages are in the system library.
  # Skip renv activation entirely — it would sandbox away the system library
  # and look for packages in an empty project library.
  if (nzchar(Sys.getenv("EPA_DOCKER"))) return()

  # Activate renv from one of the known locations:
  for (p in c("renv/activate.R",
              "/home/rstudio/renv/activate.R")) {
    if (file.exists(p)) { source(p); break }
  }
})
