#!/usr/bin/env Rscript
# Self-contained: reads original draft, modifies for PDF, renders
# Can be run via: Rscript render_pdf.R
# Or via URL:     Rscript https://raw.githubusercontent.com/neurogenomics/Empirical_Power_Analysis/main/render_pdf.R

setwd("/mnt/data/chrisc/Power_Analysis_Paper/Empirical_Power_Analysis")

message("=== Step 1: Reading original draft ===")
lines <- readLines("nature_comms_draft.rmd")

message("=== Step 2: Replacing YAML header for PDF output ===")
# Find YAML boundaries
yaml_ends <- which(lines == "---")
if (length(yaml_ends) < 2) stop("Could not find YAML front matter")

new_yaml <- c(
  "---",
  'title: "Single-Cell Study Designs Are Systematically Underpowered for Small-Effect Genes"',
  'author: "Chris Crampton, Salman Fawad, Hiru Dash, Nathan Skene"',
  'date: "`r Sys.Date()`"',
  "output:",
  "  bookdown::pdf_document2:",
  "    latex_engine: xelatex",
  "    md_extensions: +tex_math_dollars",
  "    number_sections: true",
  "    toc: false",
  "    keep_tex: true",
  "header-includes:",
  "  - \\usepackage{float}",
  "  - \\floatplacement{figure}{H}",
  "  - \\floatplacement{table}{H}",
  "editor_options:",
  "  chunk_output_type: inline",
  "bibliography: references.bib",
  "csl: https://www.zotero.org/styles/nature",
  "---"
)

# Replace YAML, keep everything after second ---
body <- lines[(yaml_ends[2] + 1):length(lines)]

# Insert latex-setup chunk right after YAML
latex_setup <- c(
  "",
  '```{r latex-setup, include=FALSE}',
  "knitr::opts_chunk$set(",
  '  out.width = "75%",',
  '  fig.align = "center",',
  '  fig.pos = "H",',
  '  out.extra = ""',
  ")",
  "```",
  ""
)

all_lines <- c(new_yaml, latex_setup, body)

# Fix out.width="99%" -> "75%" in chunk headers
all_lines <- gsub('out.width\\s*=\\s*"99%"', 'out.width="75%"', all_lines)

# Fix underscores in fig.cap strings (n_cells -> n\\_cells)
# Only in lines that contain fig.cap
fig_cap_idx <- grep("fig\\.cap", all_lines)
for (i in fig_cap_idx) {
  # Extract the fig.cap value and escape underscores
  all_lines[i] <- gsub("n_cells", "n\\\\_cells", all_lines[i])
  all_lines[i] <- gsub("n_donors", "n\\\\_donors", all_lines[i])
}

message("=== Step 3: Writing PDF version ===")
writeLines(all_lines, "nature_comms_draft_pdf.rmd")
message("Wrote nature_comms_draft_pdf.rmd")

message("=== Step 4: Rendering PDF ===")
rmarkdown::render(
  "nature_comms_draft_pdf.rmd",
  output_format = "bookdown::pdf_document2",
  output_file = "nature_comms_draft_pdf.pdf"
)
message("=== Done! PDF saved as nature_comms_draft_pdf.pdf ===")
