"""Infographic of the level rule per unit of evidence (predict_edc2.py).

Numbers measured in panel_ref.csv with the models deployed on 2026-10-27:
ERalpha has 12 deployed models that collapse into 6 units; the migration of the 83
molecule x receptor pairs comes from comparing the old rule with the new one.
"""
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch

try:
    apply_figure_style()
except NameError:
    pass
mpl.rcParams["font.family"] = "DejaVu Sans"

W, H = 100.0, 150.0
BG, INK, MUTED, RULE = "#fbfaf7", "#1d2433", "#566070", "#d5d1c7"
PHARM, SCREEN = "#1b3f73", "#6a4c93"          # pharmacological axis (ChEMBL) and screening axis (Tox21)
HIGH, MOD, LOW, NO_ALERT, NAV = "#1B3A6B", "#4F8BC9", "#A7C7E7", "#EEF3F8", "#EADCF0"
PALE = "#f1eee6"

fig = plt.figure(figsize=(11.0, 16.5), facecolor=BG)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")
TEXTS = []


# Draw a text item and register it for the overlap/overflow check.
def txt(x, y, s, size=8.5, weight="normal", color=INK, ha="left", va="center", box=None, ls=1.35):
    t = ax.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha, va=va,
                linespacing=ls, zorder=10)
    TEXTS.append((t, box))
    return t


# Rounded rectangle used as a card or panel background.
def rbox(x0, y0, x1, y1, fc, ec="none", lw=0.0, r=1.2, ls="-", z=2):
    p = FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                       boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw,
                       linestyle=ls, zorder=z)
    ax.add_patch(p)
    return p


# Short arrow between two points.
def arrow(x0, y0, x1, y1, color=None, lw=1.2, ms=7):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=ms, color=color or MUTED, lw=lw,
                                 shrinkA=0, shrinkB=0, zorder=6))


# Numbered section header with a rule underneath.
def header(y, n, title):
    ax.add_patch(Circle((4.6, y), 2.0, fc=INK, ec="none", zorder=4))
    txt(4.6, y - 0.05, str(n), size=13, weight="bold", color="white", ha="center",
        box=(2.6, y - 2.0, 6.6, y + 2.0))
    txt(8.0, y, title, size=13.5, weight="bold", box=(7.5, y - 1.7, 98.5, y + 1.7))
    ax.plot([2, 98], [y - 2.9, y - 2.9], color=RULE, lw=1.0, zorder=1)


# ------------------------------------------------------------------ title
txt(50, 145.6, "How the alert level is decided", size=20, weight="bold", ha="center",
    box=(1, 143.6, 99, 147.6))
txt(50, 141.8, "A model is not a piece of evidence. Votes are grouped before being counted.",
    size=11, color=MUTED, ha="center", box=(1, 140.4, 99, 143.2))

# ------------------------------------------------------------------ 1. collapse
header(136.5, 1, "First, repeated votes become a single piece of evidence")
txt(4, 131.6, "Real example: estrogen receptor alpha (ERα), which has 12 deployed models.",
    size=9, color=MUTED, box=(3.6, 130.6, 80, 132.6))

txt(20, 128.0, "12 MODELS", size=9.5, weight="bold", ha="center", box=(9, 127.0, 31, 129.0))
txt(74, 128.0, "6 UNITS OF EVIDENCE", size=9.5, weight="bold", ha="center",
    box=(55, 127.0, 93, 129.0))

GROUPS = [
    (PHARM, ["IC50 binding  ·  regression", "IC50 binding  ·  classification"], "ESR1_B",
     "same compounds, same features"),
    (PHARM, ["Ki binding  ·  regression", "Ki binding  ·  classification"], "ESR1_Ki_B",
     "same compounds, same features"),
    (PHARM, ["EC50 binding  ·  regression", "EC50 binding  ·  classification"], "ESR1_EC50_B",
     "same compounds, same features"),
    (PHARM, ["EC50 functional  ·  regression", "EC50 functional  ·  classification"], "ESR1_EC50_F",
     "same compounds, same features"),
    (SCREEN, ["Tox21 agonist  ·  HEK293T", "Tox21 agonist  ·  VM7"], "ERα agonist",
     "same library, different cell line"),
    (SCREEN, ["Tox21 antagonist  ·  HEK293T", "Tox21 antagonist  ·  VM7"], "ERα antagonist",
     "same library, different cell line"),
]
Y0 = 124.0
for k, (fill_color, models, unit, ratio) in enumerate(GROUPS):
    yc = Y0 - 4.1 * k
    for j, name in enumerate(models):
        yy = yc + (0.95 if j == 0 else -0.95)
        ax.add_patch(Circle((5.6, yy), 0.72, fc=fill_color, ec="white", lw=0.5, zorder=5))
        txt(7.0, yy, name, size=7.6, box=(6.4, yy - 0.95, 40, yy + 0.95))
    arrow(42.0, yc, 48.2, yc, color=fill_color)
    txt(44.5, yc + 1.6, "become 1", size=6.8, color=MUTED, ha="center",
        box=(40.5, yc + 0.8, 48.5, yc + 2.4))
    rbox(49, yc - 1.75, 76, yc + 1.75, "white", ec=fill_color, lw=1.4, r=0.9, z=3)
    ax.add_patch(Circle((51.4, yc), 0.85, fc=fill_color, ec="none", zorder=5))
    txt(53.2, yc, unit, size=8.2, weight="bold", box=(52.5, yc - 1.6, 76, yc + 1.6))
    txt(77.5, yc, ratio, size=7.2, color=MUTED, box=(77.0, yc - 1.1, 98.5, yc + 1.1))

txt(4, 97.6, "Without grouping, ERα would give 12 votes — but these are 6 pieces of evidence counted twice.",
    size=9, weight="bold", box=(3.6, 96.4, 90, 98.8))

# ------------------------------------------------------------------ 2. axes
header(92.0, 2, "Each unit belongs to one of the two axes")
for x0, fill_color, heading, body_text, cap in (
    (3, PHARM, "Pharmacological axis  (ChEMBL)",
     "Affinity and potency measured in an assay\n"
     "targeted at the receptor. One unit per\n"
     "endpoint: receptor × measure × assay kind.\n\n"
     "Positive when the CLASSIFICATION flags it.\n"
     "Regression enters as predicted potency,\n"
     "not as a separate vote.",
     "ERα 4 units  ·  ERβ 4  ·  AR 4  ·  PPARγ 5  ·  THRβ 1"),
    (51, SCREEN, "Screening axis  (Tox21)",
     "qHTS screening in cells, with an active/\n"
     "inactive label from the depositor. One unit\n"
     "per receptor and assay mode.\n\n"
     "Positive when ANY cell line flags it.\n"
     "The cell lines share up to 99.8% of the\n"
     "library, so they don't corroborate each other.",
     "ERα 2 units  ·  AR 2  ·  PPARγ 2  ·  TR 1  ·  ERβ 0"),
):
    x1 = x0 + 46
    rbox(x0, 66.0, x1, 88.0, "white", ec=fill_color, lw=1.5, r=1.2, z=2)
    rbox(x0, 83.2, x1, 88.0, fill_color, r=1.2, z=3)
    txt(x0 + 1.6, 85.6, heading, size=10, weight="bold", color="white",
        box=(x0, 83.2, x1, 88.0))
    txt(x0 + 1.6, 81.6, body_text, size=8.3, va="top", ls=1.5, box=(x0 + 0.5, 70.0, x1 - 0.5, 82.2))
    txt(x0 + 1.6, 68.6, cap, size=7.6, weight="bold", color=fill_color,
        box=(x0 + 0.5, 67.0, x1 - 0.5, 70.2))

# ------------------------------------------------------------------ 3. ladder
header(61.0, 3, "The count of units defines the level")
LEVELS = [
    (HIGH, "white", "HIGH", "the TWO positive axes", "the targeted evidence and the screening\nevidence agree — it's the strongest reading"),
    (MOD, INK, "MODERATE", "2+ positive units in a single axis", "there is weight, but it comes from a single\nexperimental origin"),
    (LOW, INK, "LOW", "exactly 1 positive unit", "isolated evidence, without corroboration"),
    (NO_ALERT, INK, "NO ALERT", "2+ units assessed, all negative", "was tested and came back negative"),
]
for k, (fc, tc, name, rule, note) in enumerate(LEVELS):
    y1 = 56.4 - 5.0 * k
    y0 = y1 - 4.5
    rbox(4, y0, 24, y1, fc, ec=RULE if fc == NO_ALERT else "none", lw=0.8, r=1.0, z=3)
    txt(14, (y0 + y1) / 2, name, size=10.5, weight="bold", color=tc, ha="center",
        box=(4, y0, 24, y1))
    txt(25.8, (y0 + y1) / 2 + 1.25, rule, size=9, weight="bold",
        box=(25.4, (y0 + y1) / 2 + 0.2, 62, (y0 + y1) / 2 + 2.3))
    txt(25.8, (y0 + y1) / 2 - 1.3, note, size=7.8, color=MUTED, ls=1.3,
        box=(25.4, y0 - 0.5, 62, (y0 + y1) / 2 + 0.1))

rbox(64, 34.6, 96.5, 56.4, NAV, ec="#b79fc6", lw=1.3, r=1.2, ls="--", z=3)
txt(80.2, 53.6, "NOT ASSESSED", size=10.5, weight="bold", color="#4A2B5A", ha="center",
    box=(65, 52.4, 95.5, 54.8))
txt(80.2, 50.8, "outside the ladder", size=8, color="#4A2B5A", ha="center",
    box=(65, 49.9, 95.5, 51.7))
txt(65.5, 48.2, "No unit was assessed, or only one\n"
                "was — and it came back negative.\n\n"
                "It is not a degree of alert: it is the absence of\n"
                "evidence. The old rule called this\n"
                "“no alert”, which is the opposite.\n\n"
                "A unit is left unassessed when its molecule\n"
                "falls outside the model's domain.",
    size=8.0, va="top", ls=1.42, color="#3b2448", box=(65, 35.2, 96.5, 48.8))

# ------------------------------------------------------------------ 4. migration
header(32.0, 4, "What changed in the 83 pairs of the reference panel")
MIG = [("high", 10, [("high", 10)]), ("moderate", 10, [("moderate", 6), ("low", 4)]),
       ("low", 23, [("low", 23)]), ("no alert", 40, [("no alert", 26), ("not assessed", 14)])]
COLORS = {"high": HIGH, "moderate": MOD, "low": LOW, "no alert": NO_ALERT, "not assessed": NAV}
TXT_COLOR = {"high": "white", "moderate": INK, "low": INK, "no alert": INK, "not assessed": "#4A2B5A"}
txt(16, 27.3, "OLD RULE", size=9, weight="bold", ha="center", box=(5, 26.4, 27, 28.2))
txt(52, 27.3, "NEW RULE", size=9, weight="bold", ha="center", box=(41, 26.4, 63, 28.2))
yy = 24.4
for name, n, destinations in MIG:
    alt = max(0.6 * n / 10 + 1.6, 2.3 * len(destinations))
    rbox(5, yy - alt, 27, yy, COLORS[name], ec=RULE if name == "no alert" else "none", lw=0.8,
         r=0.8, z=3)
    txt(16, yy - alt / 2, f"{name}  ·  {n}", size=8.6, weight="bold", color=TXT_COLOR[name], ha="center",
        box=(5, yy - alt, 27, yy))
    sub = yy
    for dn, dv in destinations:
        dalt = alt * dv / n
        arrow(27.8, sub - dalt / 2, 40.5, sub - dalt / 2, color=MUTED, lw=1.0, ms=6)
        rbox(41, sub - dalt + 0.12, 63, sub - 0.12, COLORS[dn],
             ec=RULE if dn == "no alert" else ("#b79fc6" if dn == "not assessed" else "none"),
             lw=0.8, r=0.7, z=3)
        txt(52, sub - dalt / 2, f"{dn}  ·  {dv}", size=8.2, weight="bold", color=TXT_COLOR[dn],
            ha="center", box=(41, sub - dalt - 0.3, 63, sub + 0.3))
        if dn == "not assessed":
            txt(64.5, sub - dalt / 2, "← these pairs seemed negative,\nbut no one could assess them",
                size=8, weight="bold", color="#4A2B5A", ls=1.3,
                box=(64.0, sub - dalt - 1.0, 98.5, sub + 1.0))
        sub -= dalt
    yy -= alt + 0.8

txt(4, 5.6, "The 10 “high” remained “high”, all due to convergence between the two axes.",
    size=9, weight="bold", box=(3.6, 4.4, 72, 6.8))
txt(4, 2.6, "Always read the “coverage” column next to the level: a positive fraction of 1.00 with coverage 0.17 "
            "means that only one out of seven units gave an opinion.",
    size=8.3, color=MUTED, box=(3.6, 1.4, 98, 3.8))


# After rendering, flag text that leaves its box, leaves the figure, or overlaps.
def verify():
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    out, ext = [], []
    for t, box in TEXTS:
        bb = t.get_window_extent(r)
        (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        ext.append((t, (x0, y0, x1, y1)))
        if box is not None:
            bx0, by0, bx1, by1 = box
            if x0 < bx0 - 0.1 or x1 > bx1 + 0.1 or y0 < by0 - 0.1 or y1 > by1 + 0.1:
                out.append(("outside the box", t.get_text()[:40].replace("\n", " / "),
                            tuple(round(v, 1) for v in (x0, y0, x1, y1)), box))
        if x0 < 0 or x1 > W or y0 < 0 or y1 > H:
            out.append(("outside the figure", t.get_text()[:40]))
    for i, (a, ea) in enumerate(ext):
        for b, eb in ext[i + 1:]:
            if ea[0] < eb[2] and eb[0] < ea[2] and ea[1] < eb[3] and eb[1] < ea[3]:
                out.append(("overlap", a.get_text()[:28].replace("\n", " / "),
                            b.get_text()[:28].replace("\n", " / ")))
    return out


issues = verify()
fig.savefig("alert_levels_infographic.png", dpi=220, facecolor=BG)
print(len(TEXTS), "texts |", len(issues), "issues")
for p in issues:
    print(p)
