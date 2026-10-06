# =============================================================================
# SCRIPT: compose_xgb_panels.R
# PURPOSE: Build the multi-panel XGBoost figure (fig:FigXGB in the manuscript)
#          from the five clean per-panel exports, laid out in manuscript
#          reading order and labelled A-E.
#
#          Layout (reading order A -> B, B is wide at the bottom):
#            row 1:  A  (top-left)      C  (top-right)
#            row 2:  D  (mid-left)      E  (mid-right)
#            row 3:  B  (full width, bottom)
#
#          Panel -> content -> source PNG mapping. The panel letters were
#          reassigned to match the order in which they are first referenced in
#          nature_comms_draft.rmd (panels were designed before the text was
#          written). The source PNGs keep their original descriptive names:
#
#            NEW  content                                   source PNG
#            ---  --------------------------------------    -------------------------------
#             A   Out-of-fold calibration curves            xgb_supp_calibration.png
#             B   Hyperparameter bracketing (search MAE)    xgb_supp_hp_bracketing.png
#             C   Precision/recall/F1 of class schemes      xgb_supp_prf_grid.png
#             D   Mean|SHAP| feature importance             xgb_supp_shap.png
#             E   Mean detection rate over |LFC| x donors   xgb_supp_meanpower_heatmap.png
#
#          Panel E's source PNG carries a baked-in title + subtitle that the
#          other panels do not; the top strip (TITLE_STRIP_PX rows) is cropped
#          off so the panels match.
#
# SOURCES: the five clean per-panel PNGs (no baked-in letters) in
#   Datalink/Power_Curve_Fitting/Models/xgboost/supp_panels/
#
# WRITES: the canonical composite the manuscript chunk fig:FigXGB includes:
#   Datalink/Power_Curve_Fitting/Models/xgboost/figure_xgb_composed.png
#   Re-running refreshes the draft figure in place. It does not edit the .rmd.
#
# Run (on iclneuro, where the panel PNGs and R live):
#   Rscript R/figure_composition/compose_xgb_panels.R
#   Rscript R/figure_composition/compose_xgb_panels.R --datalink /path/to/Datalink
# =============================================================================

suppressPackageStartupMessages({
  library(png)
  library(grid)
})

# ---------------------------------------------------------------------------
# Paths / args
# ---------------------------------------------------------------------------

HERE <- dirname(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  mustWork = FALSE))
if (is.na(HERE) || HERE == "") HERE <- normalizePath(".", mustWork = FALSE)
PROJECT_ROOT <- normalizePath(file.path(HERE, "..", ".."), mustWork = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[i + 1] else default
}
DATALINK <- get_arg("--datalink", file.path(PROJECT_ROOT, "Datalink"))

PANEL_DIR <- file.path(DATALINK,
  "Power_Curve_Fitting", "Models", "xgboost", "supp_panels")
OUT_FILE <- file.path(DATALINK,
  "Power_Curve_Fitting", "Models", "xgboost", "figure_xgb_composed.png")

# ---------------------------------------------------------------------------
# Layout constants
# ---------------------------------------------------------------------------

COL_W          <- 1200   # width (px) of a single-column panel
PAD            <- 110    # vertical padding below each row (px)
TITLE_STRIP_PX <- 155    # rows cropped off the top of panel E to drop its title
LETTER_SIZE    <- 66     # panel-letter font size

rd <- function(f) {
  p <- file.path(PANEL_DIR, f)
  if (!file.exists(p)) stop("Missing panel PNG: ", p)
  readPNG(p)
}

# Load panels under their NEW letters.
A <- rd("xgb_supp_calibration.png")       # A  top-left
C <- rd("xgb_supp_prf_grid.png")          # C  top-right
D <- rd("xgb_supp_shap.png")              # D  mid-left
E <- rd("xgb_supp_meanpower_heatmap.png") # E  mid-right
B <- rd("xgb_supp_hp_bracketing.png")     # B  bottom, wide

# Strip the baked-in title + subtitle off panel E.
if (nrow(E) > TITLE_STRIP_PX) {
  E <- E[(TITLE_STRIP_PX + 1):nrow(E), , , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Compose
# ---------------------------------------------------------------------------

asp <- function(im) dim(im)[1] / dim(im)[2]   # height / width
CW  <- 2 * COL_W                              # full canvas width

hA <- asp(A) * COL_W; hC <- asp(C) * COL_W
hD <- asp(D) * COL_W; hE <- asp(E) * COL_W
hB <- asp(B) * CW

h1 <- max(hA, hC) + PAD   # row 1 height
h2 <- max(hD, hE) + PAD   # row 2 height
h3 <- hB + PAD            # row 3 height
CH <- h1 + h2 + h3

dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
png(OUT_FILE, width = CW, height = CH, res = 72, bg = "white")
pushViewport(viewport(xscale = c(0, CW), yscale = c(0, CH)))

draw <- function(im, cx, top_y, w, h) {
  grid.raster(im,
    x = unit(cx, "native"), y = unit(top_y - h / 2, "native"),
    width = unit(w, "native"), height = unit(h, "native"))
}
lett <- function(s, lx, ty) {
  grid.text(s, x = unit(lx, "native"), y = unit(ty, "native"),
    just = c("left", "top"), gp = gpar(fontsize = LETTER_SIZE, fontface = "bold"))
}

# row 1: A (left), C (right)
draw(A, COL_W / 2,          CH - PAD, COL_W, hA); lett("A", 14,          CH - 8)
draw(C, COL_W + COL_W / 2,  CH - PAD, COL_W, hC); lett("C", COL_W + 14,  CH - 8)

# row 2: D (left), E (right)
t2 <- CH - h1
draw(D, COL_W / 2,          t2 - PAD, COL_W, hD); lett("D", 14,          t2 - 8)
draw(E, COL_W + COL_W / 2,  t2 - PAD, COL_W, hE); lett("E", COL_W + 14,  t2 - 8)

# row 3: B (wide)
draw(B, CW / 2,             h3 - PAD, CW,    hB); lett("B", 14,          h3 - 8)

dev.off()
cat(sprintf("WROTE %s (%dx%d)\n", OUT_FILE, CW, round(CH)))
