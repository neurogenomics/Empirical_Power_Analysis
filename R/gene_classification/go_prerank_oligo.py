# ==============================================================================
# SCRIPT: go_prerank_oligo.py
# PURPOSE: Threshold-FREE GO enrichment of oligodendrocyte gene detectability.
#
#   This is the Python half of the two-step GO pipeline (Chris, "rank in R,
#   prerank in Python"). rank_oligo_power.R ranks EVERY gene by predicted power
#   (reusing the renv-pinned xgboost + shipped predict_power.R bundle) and writes
#   the .rnk lists; this script only runs pre-ranked GSEA (`gseapy.prerank`, the
#   fgsea algorithm) over those full ranked lists. No model is loaded here.
#
#   REPLACES the previous over-representation analysis (ORA), which took the
#   top-200 easiest + bottom-200 hardest genes and ran `gseapy.enrichr` on each
#   discrete list. Nathan flagged the top/bottom-200 cutoff as arbitrary. Prerank
#   GSEA asks whether each GO term concentrates among the easy or hard end of the
#   FULL ranking. No cutoff, no discrete hit list, no background set.
#
#   Sign convention: high rank score = easy to detect. So a GO term with
#   NES > 0 is enriched among EASILY detected genes; NES < 0 is enriched among
#   HARD-to-detect genes.
#
# GENE SET (vendored offline for reproducibility):
#   Datalink/Gene_Classification/gene_sets/GO_Biological_Process_2023.gmt
#   The .gmt is committed to the repo so the container needs NO runtime network
#   (gseapy.prerank never downloads from Enrichr). Results depend on a pinned
#   snapshot, not on whatever Enrichr serves that day.
#
# TWO GENE SETS (both requested), read from rank_oligo_power.R:
#   allgene : the full oligo transcriptome.
#   pc      : protein-coding only.
#
# INPUTS (Datalink/Gene_Classification/GO_Enrichment/, from rank_oligo_power.R):
#   oligo_<set>_prerank_input.rnk        symbol<TAB>score, no header, descending
#
# OUTPUTS (Datalink/Gene_Classification/GO_Enrichment/):
#   oligo_<set>_go_prerank_full.csv      every tested GO term
#   oligo_<set>_go_prerank_sig.csv       terms passing FDR_Q_SIG
#   oligo_<set>_go_prerank_nes.png/.pdf  top easy vs hard terms by NES
#       where <set> in {allgene, pc}
#
# RUN (after rank_oligo_power.R; needs gseapy + matplotlib, both pinned in the
# Docker image's python-requirements.txt):
#   python3 R/gene_classification/go_prerank_oligo.py
# ==============================================================================

import os
import sys

import pandas as pd

# ------------------------------------------------------------------ PARAMETERS
ROOT = os.environ.get(
    "POWER_ANALYSIS_ROOT",
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")),
)

GO_DIR = os.path.join(ROOT, "Datalink", "Gene_Classification", "GO_Enrichment")
GMT_PATH = os.path.join(
    ROOT, "Datalink", "Gene_Classification", "gene_sets",
    "GO_Biological_Process_2023.gmt",
)

SETS = ["allgene", "pc"]          # both requested; oligo only

# GSEA settings.
MIN_GENE_SET   = 15
MAX_GENE_SET   = 500
PERMUTATIONS   = 1000
SEED           = 42
THREADS        = 4
FDR_Q_SIG      = 0.25          # "sig" file cutoff (standard GSEA FDR)
N_PLOT_PER_DIR = 15            # top-N easy and top-N hard terms in the NES plot


# --------------------------------------------------------------------- HELPERS
def log(msg):
    print(f"[go_prerank] {msg}", flush=True)


def load_rnk(tag):
    path = os.path.join(GO_DIR, f"oligo_{tag}_prerank_input.rnk")
    if not os.path.isfile(path):
        sys.exit(f"missing ranked input '{path}' -- run rank_oligo_power.R first")
    rnk = pd.read_csv(path, sep="\t", header=None, names=["symbol", "score"])
    rnk = rnk.dropna(subset=["symbol"])
    rnk["score"] = pd.to_numeric(rnk["score"], errors="coerce")
    rnk = rnk.dropna(subset=["score"]).reset_index(drop=True)
    return rnk


def run_prerank(rnk, tag):
    import gseapy as gp
    log(f"[{tag}] prerank on {len(rnk)} genes ...")
    pre = gp.prerank(
        rnk=rnk, gene_sets=GMT_PATH,
        min_size=MIN_GENE_SET, max_size=MAX_GENE_SET,
        permutation_num=PERMUTATIONS, seed=SEED, threads=THREADS,
        outdir=None, no_plot=True, verbose=True,
    )
    res = pre.res2d.copy()
    # Normalise the FDR/NES column names across gseapy versions.
    ren = {}
    for c in res.columns:
        cl = c.lower()
        if cl in ("fdr q-val", "fdr", "fdr_q-val"):
            ren[c] = "FDR"
        elif cl in ("nom p-val", "nom_p-val", "p-val", "pval"):
            ren[c] = "NOM_p"
        elif cl == "nes":
            ren[c] = "NES"
        elif cl == "es":
            ren[c] = "ES"
    res = res.rename(columns=ren)
    for col in ("NES", "FDR", "NOM_p", "ES"):
        if col in res.columns:
            res[col] = pd.to_numeric(res[col], errors="coerce")
    res = res.sort_values("NES", ascending=False).reset_index(drop=True)
    return res


def plot_nes(res, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sig = res[res["FDR"] < FDR_Q_SIG].copy() if "FDR" in res.columns else res.copy()
    if sig.empty:
        log(f"[{tag}] no terms below FDR {FDR_Q_SIG}; skipping NES plot")
        return
    easy = sig[sig["NES"] > 0].nlargest(N_PLOT_PER_DIR, "NES")
    hard = sig[sig["NES"] < 0].nsmallest(N_PLOT_PER_DIR, "NES")
    sel = pd.concat([hard, easy]).sort_values("NES")
    if sel.empty:
        return
    term_col = "Term" if "Term" in sel.columns else sel.columns[0]
    labels = sel[term_col].astype(str).str.slice(0, 60)
    colours = ["#0072B2" if v > 0 else "#D55E00" for v in sel["NES"]]
    fig, ax = plt.subplots(figsize=(9, 0.34 * len(sel) + 1.2))
    ax.barh(range(len(sel)), sel["NES"], color=colours)
    ax.set_yticks(range(len(sel)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("NES  (>0 easier to detect,  <0 harder to detect)")
    ax.set_title(f"Oligo {tag}: GO terms by detectability (prerank GSEA, FDR<{FDR_Q_SIG})",
                 fontsize=10)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(GO_DIR, f"oligo_{tag}_go_prerank_nes.{ext}"), dpi=200)
    plt.close(fig)
    log(f"[{tag}] wrote NES plot ({len(sel)} terms)")


# ------------------------------------------------------------------------ MAIN
def main():
    if not os.path.isfile(GMT_PATH):
        sys.exit(f"missing vendored gene set '{GMT_PATH}'")
    log(f"gene set: {GMT_PATH}")

    for tag in SETS:
        rnk = load_rnk(tag)
        log(f"[{tag}] {len(rnk)} ranked symbols "
            f"(score range {rnk['score'].min():.4f}..{rnk['score'].max():.4f})")

        res = run_prerank(rnk, tag)
        res.to_csv(os.path.join(GO_DIR, f"oligo_{tag}_go_prerank_full.csv"), index=False)
        if "FDR" in res.columns:
            sig = res[res["FDR"] < FDR_Q_SIG]
            sig.to_csv(os.path.join(GO_DIR, f"oligo_{tag}_go_prerank_sig.csv"), index=False)
            log(f"[{tag}] {len(res)} terms tested; {len(sig)} at FDR<{FDR_Q_SIG}")
        plot_nes(res, tag)

    log("done.")


if __name__ == "__main__":
    main()
