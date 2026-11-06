"""Figure of the fair test: paired differences between training conditions.

Panels: RMSE gain (regression), ROC-AUC gain on ChEMBL and on Tox21 (classification), for conditions B and C against A,
on the same test set. Reads the long table written by analyses/fair_test_scaffold/fair_test_analysis.py.

Usage: python build_fig_fair_test_differences.py <fair_long.csv> [output.png]
"""
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# Style helpers: apply_figure_style and panel_letter.
META_GREY = "#888888"


# Apply the shared publication style (open frame, 8/7/6 pt sizes) to matplotlib rcParams.
def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


# Write a bold lowercase panel letter (a, b, c) above the top-left corner of an axes.
def panel_letter(ax, letter, dx=-0.18, dy=1.02, case="lower", fontsize=None):
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    ax.text(dx, dy, s, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")



apply_figure_style(sizes=(8, 7, 6))

CSV = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else "fig_fair_test_differences.png"
D = pd.read_csv(CSV)
D = D[~((D.group != "tox21") & (D.random_seed == 15))]

C = ["A_no_large_groups", "B_with_large_groups", "C_equal_size"]

# Exclude numerical failures (R2 < -5 in any condition for the same endpoint and seed)
reg = D[D.task == "regression"].copy()
failure = reg.groupby(["group", "endpoint", "random_seed"]).r2.transform("min") < -5
regok = reg[~failure]

cl = D[D.task == "classification"]

BLUE = "#1f6fb4"
ORG = "#d9822b"


# Per-endpoint median difference of condition B and C against A, signed so that positive means better than A.
def deltas(d, metric, groups, lower_better):
    m = d[d.group.isin(groups)].groupby(["endpoint", "condition"])[metric].median().unstack().dropna(subset=C)
    s = -1 if lower_better else 1       # positive = better than A
    return {"B": (s * (m["B_with_large_groups"] - m["A_no_large_groups"])).values,
            "C": (s * (m["C_equal_size"] - m["A_no_large_groups"])).values}


panels = [("Regression\nRMSE gain", deltas(regok, "rmse", ["ic50", "layers"], True), "RMSE gain (pAct)"),
          ("ChEMBL classification\nROC-AUC gain", deltas(cl, "roc_auc", ["ic50", "layers"], False), "ROC-AUC gain"),
          ("Tox21 classification\nROC-AUC gain", deltas(cl, "roc_auc", ["tox21"], False), "ROC-AUC gain")]

fig, axs = plt.subplots(1, 3, figsize=(7.4, 3.3), gridspec_kw=dict(wspace=0.62))
rg = np.random.default_rng(1)
for ax, (tit, dd, yl) in zip(axs, panels):
    for k, (key, col) in enumerate((("B", BLUE), ("C", ORG))):
        v = dd[key]
        x = np.full(len(v), k) + rg.uniform(-0.13, 0.13, len(v))
        ax.scatter(x, v, s=11, color=col, edgecolor="white", linewidth=0.4, alpha=0.9, zorder=3)
        bs = [np.median(rg.choice(v, len(v), replace=True)) for _ in range(4000)]
        lo, hi = np.percentile(bs, [2.5, 97.5])
        med = np.median(v)
        ax.plot([k - 0.25, k + 0.25], [med, med], color="black", lw=1.6, zorder=4, solid_capstyle="butt")
        ax.plot([k, k], [lo, hi], color="black", lw=1.0, zorder=4)
    ax.axhline(0, color="black", lw=0.6, ls=(0, (3, 2)), zorder=0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["B", "C"])
    ax.set_xlim(-0.6, 1.6)
    ax.set_title(tit, loc="left")
    ax.set_ylabel(yl)
    ax.set_xlabel(f"n = {len(dd['B'])} endpoints")
for a, l in zip(axs, "abc"):
    panel_letter(a, l, dx=-0.30)
fig.text(0.5, -0.10,
         "Gain > 0 = better than A, on the same test set.  A = large groups out of training (60% of the data)  B = large groups in training (80%)\n"
         "C = B subsampled to the size of A.  Black bar = median; vertical line = 95% bootstrap CI across endpoints.",
         ha="center", va="top", fontsize=6, color="#444444")
fig.savefig(OUT, dpi=300, bbox_inches="tight")
for tit, dd, _ in panels:
    print(tit.replace("\n", " "), "| n =", len(dd["B"]), "| median B = %.4f, C = %.4f" % (np.median(dd["B"]), np.median(dd["C"])))
