"""Figure for the BindingDB layer merging experiment (literature).

One sentence: adding the literature fraction does not change performance in the
ChEMBL chemical space and transforms it in the new chemical space.

Data: experiment_merge.csv, 8 endpoints x 3 scaffold-split seeds, with the
identical test in both conditions and fixed algorithm (RandomForest, 400 trees).
Usage: python build_fig_merge_bindingdb.py <experiment_merge.csv>
"""
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# --- style: three sizes mapped to paper (base / annotation / tick)
BASE, ANNOT, TICK = 8, 7, 6
mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": BASE,
    "axes.titlesize": BASE, "axes.labelsize": BASE,
    "xtick.labelsize": TICK, "ytick.labelsize": TICK,
    "legend.fontsize": ANNOT, "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "savefig.facecolor": "white", "figure.facecolor": "white",
})
CURRENT, NEW, INK, MUTED, RULE = "#8c8c8c", "#1b3f73", "#1d2433", "#5a6472", "#d9d6cf"

PATH = sys.argv[1] if len(sys.argv) > 1 else "experiment_merge.csv"
T = pd.read_csv(PATH)

METRICS = [("roc_auc", "ROC-AUC", 0.5, "random"),
            ("pr_auc", "PR-AUC", None, None),
            ("r2", "R²", 0.0, "predicts the average")]
ROWS = [("ChEMBL_test", "Assessed in the ChEMBL chemical space"),
          ("BDB_test", "Assessed in the new chemical space, from BindingDB")]

# median of the 3 seeds per endpoint, for each condition
med = (T.groupby(["endpoint", "test", "condition"])[[m for m, *_ in METRICS]]
         .median().reset_index())

fig, axes = plt.subplots(2, 3, figsize=(6.9, 4.6))
fig.subplots_adjust(left=0.088, right=0.975, top=0.790, bottom=0.135,
                    hspace=0.40, wspace=0.40)
TEXTS = []
# R2 axis floor: one endpoint reaches -3.97 and would compress the entire panel,
# so it appears as an out-of-scale marker at the floor, with the value written next to it
R2_FLOOR = -1.25


# Register a text artist so the geometric check can walk every string.
def txt(obj):
    TEXTS.append(obj)
    return obj


for i, (test, header) in enumerate(ROWS):
    sub = med[med.test == test]
    for j, (col, label, zero, zero_txt) in enumerate(METRICS):
        ax = axes[i, j]
        d = sub.pivot_table(index="endpoint", columns="condition", values=col).dropna()
        floor = R2_FLOOR if col == "r2" and d.min().min() < R2_FLOOR else None
        out_of_scale = []
        # a thin line per endpoint
        for ep, r in d.iterrows():
            a, b = r.A, r.B
            if floor is not None and (a < floor or b < floor):
                ax.plot([0, 1], [max(a, floor), max(b, floor)], color=MUTED, lw=0.7,
                        alpha=0.55, zorder=2)
                for x, v in ((0, a), (1, b)):
                    if v < floor:
                        out_of_scale.append(v)
                        # shifts each triangle so that two out-of-scale points on the
                        # same x don't fall exactly on top of each other
                        ax.plot([x + 0.05 * (len(out_of_scale) - 1)], [floor], marker="v",
                                ms=4.2, color=MUTED, mec="white", mew=0.5, zorder=3)
                    else:
                        ax.plot([x], [v], marker="o", ms=2.6, color=MUTED, alpha=0.55,
                                mfc="white", mew=0.6, zorder=2)
                continue
            ax.plot([0, 1], [a, b], color=MUTED, lw=0.7, alpha=0.55,
                    marker="o", ms=2.6, mfc="white", mew=0.6, zorder=2)
        # median across endpoints, highlighted
        ax.plot([0, 1], [d.A.median(), d.B.median()], color=NEW, lw=2.0,
                marker="o", ms=5.0, mfc=NEW, mec="white", mew=0.8, zorder=4,
                solid_capstyle="round")
        if zero is not None:
            ax.axhline(zero, color=INK, lw=0.6, ls=(0, (3, 2)), zorder=1)
            txt(ax.text(0.97, zero, zero_txt, size=TICK, color=MUTED, va="top",
                        ha="right", transform=ax.get_yaxis_transform(), zorder=5))
        ax.set_xlim(-0.28, 1.28)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["current", "+ BindingDB"] if i == 1 else [])
        if i == 1:
            for lab, fill_color in zip(ax.get_xticklabels(), (CURRENT, NEW)):
                lab.set_color(fill_color)
                lab.set_fontweight("bold")
        txt(ax.set_ylabel(label, labelpad=2))
        ax.margins(y=0.16)
        if floor is not None:
            ax.set_ylim(bottom=floor - 0.18)
            if out_of_scale:
                band = f"{min(out_of_scale):.1f} to {max(out_of_scale):.1f}"
                txt(ax.text(0.5, 0.012, f"▾ {len(out_of_scale)} outside scale ({band})",
                            transform=ax.transAxes, ha="center", va="bottom",
                            size=TICK, color=MUTED)) if len(out_of_scale) > 1 else \
                txt(ax.text(0.5, 0.012, f"▾ outside scale ({band.split(' to ')[0]})",
                            transform=ax.transAxes, ha="center", va="bottom",
                            size=TICK, color=MUTED))
        ax.grid(axis="y", color=RULE, lw=0.45, zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        # median variation, the only result annotation in the panel
        diff = d.B.median() - d.A.median()
        txt(ax.annotate(f"{diff:+.2f}", xy=(0.5, 0.965),
                        xycoords="axes fraction", ha="center", va="top",
                        size=ANNOT, color=NEW if diff > 0 else CURRENT, fontweight="bold"))
    # band header, anchored at the top of the already positioned panel
    txt(fig.text(0.012, axes[i, 0].get_position().y1 + 0.016, header, size=BASE,
                 fontweight="bold", color=INK, ha="left", va="bottom"))

txt(fig.text(0.012, 0.995, "Add the BindingDB literature fraction to training",
             size=BASE + 2.5, fontweight="bold", color=INK, ha="left", va="top"))
txt(fig.text(0.012, 0.950,
             "Same test set in both conditions; fixed algorithm and hyperparameters. "
             "Median of 8 endpoints × 3 seeds.",
             size=ANNOT, color=MUTED, ha="left", va="top"))
txt(fig.text(0.975, 0.995, "higher = better", size=ANNOT, color=MUTED,
             ha="right", va="top", style="italic"))

fig.legend(handles=[Line2D([], [], color=MUTED, lw=0.7, marker="o", ms=2.6, mfc="white",
                           label="one endpoint"),
                    Line2D([], [], color=NEW, lw=2.0, marker="o", ms=5.0,
                           label="median of the 8 endpoints")],
           loc="lower center", bbox_to_anchor=(0.54, 0.002), ncol=2,
           frameon=False, handletextpad=0.5, columnspacing=1.6)


# Geometric check: no visible text box may overlap another or leave the figure.
def verify():
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    # ticks outside the axis limits still exist as Text and do not appear in the
    # render; counting them would produce a phantom overlap with whatever is above the panel
    inside = set()
    for ax_ in fig.axes:
        y0, y1 = sorted(ax_.get_ylim())
        inside |= {id(t) for t in ax_.get_yticklabels()
                   if not (y0 <= t.get_position()[1] <= y1)}
    boxes = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
              if t.get_text().strip() and t.get_visible() and id(t) not in inside]
    probs = [(a.get_text()[:26], b.get_text()[:26])
             for i, (a, ba) in enumerate(boxes) for b, bb in boxes[i + 1:]
             if ba.overlaps(bb)]
    outside = [t.get_text()[:26] for t, bb in boxes
            if not (fig.bbox.x0 <= bb.x0 and bb.x1 <= fig.bbox.x1
                    and fig.bbox.y0 <= bb.y0 and bb.y1 <= fig.bbox.y1)]
    return probs, outside


probs, outside = verify()
fig.savefig("fig_merge_bindingdb.png", dpi=300, bbox_inches="tight")
print(f"{len(boxes) if (boxes := None) else len(TEXTS)} texts registered | "
      f"overlaps: {len(probs)} | outside the figure: {len(outside)}")
for p in probs[:8]:
    print("  overlaps:", p)
for f in outside[:8]:
    print("  outside:", f)
