# =============================================================================
# SCRIPT: compose_xgb_panels.R
# PURPOSE: Build the multi-panel XGBoost figure (fig:FigXGB in the manuscript)
#          from the five clean per-panel exports, laid out in manuscript
#          reading order and labelled A-E.
#
#          Layout (reading order A -> E):
#            left column:   A (tall), then D directly below it
#            right column:  B, then C, then E
#          The two columns stack independently, so D rises to sit under A
#          instead of waiting for the taller right column (B + C) to end.
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
# Run (on a host where the panel PNGs and R live):
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
B <- rd("xgb_supp_hp_bracketing.png")     # B  top-right, above C

# Strip the baked-in title + subtitle off panel E.
if (nrow(E) > TITLE_STRIP_PX) {
  E <- E[(TITLE_STRIP_PX + 1):nrow(E), , , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Compose
# ---------------------------------------------------------------------------

asp <- function(im) dim(im)[1] / dim(im)[2]   # height / width
CW  <- 2 * COL_W                              # full canvas width

hA <- asp(A) * COL_W; hB <- asp(B) * COL_W; hC <- asp(C) * COL_W
hD <- asp(D) * COL_W; hE <- asp(E) * COL_W

# Each column stacks independently, so D sits directly under A and E
# directly under C (no shared row 2 boundary leaving a gap under the
# shorter column). Positions below are measured down from the top edge.
topA <- PAD;                 topD <- topA + hA + PAD
topB <- PAD;                 topC <- topB + hB + PAD
topE <- topC + hC + PAD
CH   <- max(topD + hD, topE + hE)

dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
png(OUT_FILE, width = CW, height = CH, res = 72, bg = "white")
pushViewport(viewport(xscale = c(0, CW), yscale = c(0, CH)))

draw <- function(im, cx, top, w, h) {
  grid.raster(im,
    x = unit(cx, "native"), y = unit(CH - top - h / 2, "native"),
    width = unit(w, "native"), height = unit(h, "native"))
}
lett <- function(s, lx, top) {
  grid.text(s, x = unit(lx, "native"), y = unit(CH - top, "native"),
    just = c("left", "top"), gp = gpar(fontsize = LETTER_SIZE, fontface = "bold"))
}

# Each letter sits PAD - 8 px above the top of its panel.
lx_l <- 14; lx_r <- COL_W + 14; cx_l <- COL_W / 2; cx_r <- COL_W + COL_W / 2
draw(A, cx_l, topA, COL_W, hA); lett("A", lx_l, topA - PAD + 8)
draw(D, cx_l, topD, COL_W, hD); lett("D", lx_l, topD - PAD + 8)
draw(B, cx_r, topB, COL_W, hB); lett("B", lx_r, topB - PAD + 8)
draw(C, cx_r, topC, COL_W, hC); lett("C", lx_r, topC - PAD + 8)
draw(E, cx_r, topE, COL_W, hE); lett("E", lx_r, topE - PAD + 8)

dev.off()
cat(sprintf("WROTE %s (%dx%d)\n", OUT_FILE, CW, round(CH)))
