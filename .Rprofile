local({
  # 1. If a Docker image has installed packages outside the repo at
  #    RENV_PATHS_LIBRARY, prepend it to .libPaths() so we always see the
  #    image's pinned library even when this repo is bind-mounted on top.
  #    Strictly additive: a no-op when RENV_PATHS_LIBRARY isn't set, so this
  #    is safe to roll out before the new Docker image is in production.
  image_lib <- Sys.getenv("RENV_PATHS_LIBRARY", unset = "")
  if (nzchar(image_lib) && dir.exists(image_lib)) {
    cur <- .libPaths()
    if (!(image_lib %in% cur)) {
      .libPaths(c(image_lib, cur))
    }
  }

  # 2. Activate renv from one of the known locations:
  #      - "renv/activate.R"               (normal host workflow, repo CWD)
  #      - "/home/rstudio/renv/activate.R" (legacy: pre-rewrite Docker image
  #                                          baked it under rstudio's home)
  #      - "/opt/build/renv/activate.R"    (post-rewrite Docker image baked it
  #                                          to /opt/build during install)
  for (p in c("renv/activate.R",
              "/home/rstudio/renv/activate.R",
              "/opt/build/renv/activate.R")) {
    if (file.exists(p)) { source(p); break }
  }
})
