"""Infographic of the endocrine disruptor predictor (multi-layer consensus).

Numbers verified in the user's pipeline:
  deployed models: 18 (IC50) + 32 (Ki 11, Kd 2, EC50 19) + 11 (Tox21) = 61  (count of the .joblib files)
  voting and alert rules: predict_edc2.py; external BindingDB: 817 molecule-target pairs, 769 molecules.
"""
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch, RegularPolygon

try:
    apply_figure_style()
except NameError:
    pass
mpl.rcParams["font.family"] = "DejaVu Sans"

W, H = 100.0, 144.0
FIGW, FIGH = 12.0, 17.28            # 0.12 in per unit on both axes
BG, INK, MUTED, RULE = "#fbfaf7", "#1d2433", "#566070", "#d5d1c7"
C_IC50, C_KEC, C_TOX, C_BDB = "#1b3f73", "#2b7bb9", "#6a4c93", "#6b7280"
HIGH, MOD, BAI, NO_ALERT = "#b3261e", "#ee8a3a", "#f3c75a", "#cfd6e0"
NEUTRAL = "#e6eaf0"

fig = plt.figure(figsize=(FIGW, FIGH), facecolor=BG)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")

TEXTS = []


# Draw a text item and register it for the overlap/overflow check.
def txt(x, y, s, size=8.5, weight="normal", color=INK, ha="left", va="center", style="normal", box=None, ls=1.3):
    t = ax.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha, va=va, fontstyle=style,
                linespacing=ls, zorder=10)
    TEXTS.append((t, box))
    return t


# Rounded rectangle used as a card or panel background.
def rbox(x0, y0, x1, y1, fc, ec="none", lw=0.0, r=1.2, ls="-", z=2):
    p = FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw, ls=ls, zorder=z)
    ax.add_patch(p)
    return p


# Arrow between two points in figure coordinates.
def arrow(x0, y0, x1, y1, color=MUTED, lw=2.0, ms=14):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=ms, lw=lw, color=color,
                                 zorder=6, shrinkA=0, shrinkB=0))


# Rounded pill with a centred bold label.
def chip(x0, y0, w, h, label, fc, tc="white", size=7.5, ec="none", ls="-"):
    rbox(x0, y0, x0 + w, y0 + h, fc, ec=ec, lw=1.1 if ec != "none" else 0, r=h / 2, ls=ls, z=4)
    txt(x0 + w / 2, y0 + h / 2, label, size=size, weight="bold", color=tc, ha="center", box=(x0, y0, x0 + w, y0 + h))


# Numbered section header with a horizontal rule underneath.
def header(y, n, title):
    ax.add_patch(Circle((4.6, y), 2.0, fc=INK, ec="none", zorder=4))
    txt(4.6, y - 0.05, str(n), size=13, weight="bold", color="white", ha="center", box=(2.6, y - 2.0, 6.6, y + 2.0))
    txt(8.0, y, title, size=14, weight="bold", box=(7.5, y - 1.7, 98.5, y + 1.7))
    ax.plot([2, 98], [y - 2.8, y - 2.8], color=RULE, lw=1.0, zorder=1)


# Small molecule glyph standing for the input compound.
def mol_icon(cx, cy, r=2.0):
    ax.add_patch(RegularPolygon((cx, cy), 6, radius=r, orientation=0, fc="none", ec=INK, lw=1.8, zorder=5))
    ax.add_patch(Circle((cx, cy), r * 0.55, fc="none", ec=INK, lw=1.0, zorder=5))
    for ang, d, col in ((30, (1.9, 1.1), "#c0392b"), (210, (-1.9, -1.1), "#2c6fbb")):
        a = np.deg2rad(ang)
        v = (cx + r * np.cos(a), cy + r * np.sin(a))
        e = (v[0] + d[0], v[1] + d[1])
        ax.plot([v[0], e[0]], [v[1], e[1]], color=INK, lw=1.6, zorder=5, solid_capstyle="round")
        ax.add_patch(Circle(e, 0.75, fc=col, ec="white", lw=0.5, zorder=6))
    ax.plot([cx, cx], [cy - r, cy - r - 1.7], color=INK, lw=1.6, zorder=5, solid_capstyle="round")


# Small check mark used in the quality-control list.
def tick(cx, cy):
    ax.add_patch(Circle((cx, cy), 0.95, fc=INK, ec="none", zorder=4))
    ax.plot([cx - 0.45, cx - 0.1, cx + 0.5], [cy, cy - 0.4, cy + 0.45], color="white", lw=1.4, zorder=5,
            solid_capstyle="round")


# ---------------------------------------------------------------- title
txt(50, 140.2, "How the model identifies possible endocrine disruptors", size=20, weight="bold",
    ha="center", box=(1, 138.3, 99, 142.1))
txt(50, 136.3, "A molecule enters (SMILES). Up to 61 models vote. A single alert level comes out per receptor.",
    size=11, color=MUTED, ha="center", box=(1, 134.9, 99, 137.7))

# ---------------------------------------------------------------- 1. data
header(131.2, 1, "Data: where they come from and what they're for")
cards = [
    (C_IC50, "ChEMBL 37  ·  IC50",
     "Potency of the molecule at the receptor\n(IC50), in human protein.\nBinding (B) and functional (F)\nare in separate endpoints.",
     "9 endpoints  ·  18 models", True, False),
    (C_KEC, "ChEMBL 37  ·  Ki, Kd, EC50",
     "Affinity (Ki, Kd) and effect on\ncell (EC50). Each measure has\nits own dataset, also separated\ninto binding and functional.",
     "17 endpoints  ·  32 models", True, False),
    (C_TOX, "Tox21  ·  AC50",
     "qHTS screening across 12 assays.\nLabel: active or inactive.\nExcluded are actives that may be\ncytotoxicity (antagonists)\nor autofluorescence (BLA assays).",
     "12 assays  ·  11 models", True, False),
    (C_BDB, "BindingDB",
     "IC50 records for molecules\nnot in ChEMBL: 817\nmolecule-target pairs (769 molecules).\nTests the IC50 models.",
     "does not train  ·  does not vote", False, True),
]
CY1, CY0 = 126.6, 106.4
for i, (col, title, body, count, voting, dashed) in enumerate(cards):
    x0 = 2 + i * 24.5
    x1 = x0 + 22.5
    rbox(x0, CY0, x1, CY1, fc="white", ec=col, lw=1.5, r=1.2, ls="--" if dashed else "-", z=2)
    strip = FancyBboxPatch((x0, CY1 - 5.1), 22.5, 5.1, boxstyle="round,pad=0,rounding_size=1.2", fc=col,
                           ec="none", zorder=3)
    ax.add_patch(strip)
    strip.set_clip_path(Rectangle((x0, CY1 - 3.9), 22.5, 3.9, transform=ax.transData))
    txt(x0 + 1.2, CY1 - 1.95, title, size=9.5, weight="bold", color="white", box=(x0, CY1 - 3.9, x1, CY1))
    txt(x0 + 1.2, CY1 - 5.5, body, size=8.5, va="top", ls=1.55, box=(x0 + 0.4, CY0 + 6.0, x1 - 0.4, CY1 - 3.9))
    txt(x0 + 1.2, CY0 + 5.4, count, size=8.5, weight="bold", color=col, box=(x0 + 0.4, CY0 + 4.0, x1 - 0.4, CY0 + 6.8))
    if voting:
        chip(x0 + 1.2, CY0 + 1.2, 6.2, 2.4, "TRAINS", INK)
        chip(x0 + 8.2, CY0 + 1.2, 11.6, 2.4, "VOTES IN CONSENSUS", INK)
    else:
        chip(x0 + 1.2, CY0 + 1.2, 9.4, 2.4, "VALID ONLY", "white", tc=C_BDB, ec=C_BDB, ls="--")
txt(2, 104.1, "The layers never mix: each one trains its own models. Only the votes from these models enter the consensus, not the raw data.",
    size=8.5, color=MUTED, box=(1.5, 103.2, 98.5, 105.0))
txt(2, 102.3, "PubChem BioAssay was checked and is not included: its IC50 records for these targets are copies of ChEMBL and BindingDB.",
    size=8.5, color=MUTED, box=(1.5, 101.4, 98.5, 103.2))

# ---------------------------------------------------------------- 2. cleaning
header(97.6, 2, "Cleaning: only reliable records enter")
txt(2, 92.7, "ChEMBL 37  ·  IC50, Ki, Kd, EC50", size=9.5, weight="bold", color=C_IC50, box=(1.5, 91.9, 60, 93.5))
txt(65, 92.7, "Tox21  ·  AC50", size=9.5, weight="bold", color=C_TOX, box=(64.5, 91.9, 98.5, 93.5))
steps_c = ["Single human protein, target confidence ≥ 8",
           "Exact value (=), converted to nM, plausible band",
           "Standardized structure: no salts or solvents, 100–1000 Da",
           "Median per molecule; discards replicates with more than 1 log of difference",
           "Potency becomes pAct (higher = more potent); minimum of 50 records per endpoint"]
tints_c = ["#e1ebf7", "#d0dff1", "#bfd3ec", "#aec7e6", "#9dbbe0"]
steps_t = ["12 qHTS assays; inconclusive results are excluded",
           "Cytotoxicity (antagonists) or\nautofluorescence (BLA): no label",
           "Replicates resolved by majority vote",
           "Label: active or inactive; 11 of the 12 assays\nbecome a model (TR agonist has only 27 actives)"]
tints_t = ["#ece6f3", "#e0d7ec", "#d4c8e5", "#c8badd"]
for i, (s, tcol) in enumerate(zip(steps_c, tints_c)):
    w = 58 - 2 * i
    y1 = 90.9 - 3.5 * i
    rbox(31 - w / 2, y1 - 3.0, 31 + w / 2, y1, tcol, r=0.9, z=3)
    txt(31, y1 - 1.5, s, size=8.5, ha="center", box=(31 - w / 2, y1 - 3.0, 31 + w / 2, y1))
for i, (s, tcol) in enumerate(zip(steps_t, tints_t)):
    w = 32 - 1 * i
    y1 = 90.9 - 3.5 * i
    rbox(81 - w / 2, y1 - 3.0, 81 + w / 2, y1, tcol, r=0.9, z=3)
    txt(81, y1 - 1.5, s, size=8.0 if "\n" in s else 8.5, ls=1.15, ha="center", box=(81 - w / 2, y1 - 3.0, 81 + w / 2, y1))

# ---------------------------------------------------------------- 3. training
header(70.4, 3, "Training: each endpoint (receptor × assay kind) gets two models")
T3 = 64.8
txt(9.5, T3, "Molecule", size=9.5, weight="bold", ha="center", box=(2, T3 - 1, 17, T3 + 1))
mol_icon(9.5, T3 - 6.6, r=2.3)
txt(9.5, T3 - 14.0, "SMILES\nstandardized", size=8.5, ha="center", box=(2, T3 - 16, 17, T3 - 12))
arrow(17.4, T3 - 8.5, 20.6, T3 - 8.5)
txt(32.5, T3, "Numerical representation", size=9.5, weight="bold", ha="center", box=(21, T3 - 1, 44, T3 + 1))
rbox(21, T3 - 7.2, 44, T3 - 2.8, "#dfe7f2", ec="#9aa9c0", lw=1, r=0.9, z=3)
txt(32.5, T3 - 5.0, "2,048 counts of\nsubstructures (ECFP4)", size=8.5, ha="center", box=(21, T3 - 7.2, 44, T3 - 2.8))
rbox(21, T3 - 12.6, 44, T3 - 8.2, "#ece8df", ec="#b9b2a2", lw=1, r=0.9, z=3)
txt(32.5, T3 - 10.4, "132 2D RDKit descriptors\n(mass, polarity, rings)", size=8.5, ha="center",
    box=(21, T3 - 12.6, 44, T3 - 8.2))
chip(23.0, T3 - 17.0, 19.0, 3.0, "= 2,180 numbers per molecule", INK, size=8)
arrow(44.6, T3 - 8.5, 46.4, T3 - 8.5)
txt(61.5, T3, "Two models per endpoint", size=9.5, weight="bold", ha="center", box=(47, T3 - 1, 76, T3 + 1))
rbox(47, T3 - 6.0, 76, T3 - 1.6, NEUTRAL, ec="#9aa9c0", lw=1, r=0.9, z=3)
txt(48.2, T3 - 3.0, "Regression", size=9, weight="bold", box=(47, T3 - 6.0, 76, T3 - 1.6))
txt(48.2, T3 - 4.7, "predicts the pAct value (interaction strength)", size=8.5, box=(47, T3 - 6.0, 76, T3 - 1.6))
rbox(47, T3 - 11.4, 76, T3 - 7.0, NEUTRAL, ec="#9aa9c0", lw=1, r=0.9, z=3)
txt(48.2, T3 - 8.4, "Classification", size=9, weight="bold", box=(47, T3 - 11.4, 76, T3 - 7.0))
txt(48.2, T3 - 10.1, "active if pAct ≥ 6 (potency up to 1 µM)", size=8.5, box=(47, T3 - 11.4, 76, T3 - 7.0))
txt(47.5, T3 - 12.9,
    "Algorithm selected by cross-validation\n(Random Forest, XGBoost, Ridge/SVR, Logistic/SVC).\nTox21: classification only.",
    size=8, color=MUTED, va="top", ls=1.3, box=(47, T3 - 18.0, 77, T3 - 12.0))
arrow(76.5, T3 - 8.5, 78.4, T3 - 8.5)
txt(79.5, T3, "Quality controls", size=9.5, weight="bold", box=(79, T3 - 1, 98.5, T3 + 1))
checks = ["Test on unseen\nscaffolds (80/20)", "y-scrambling: shuffled\nlabels should fail",
          "Applicability domain\n(Tanimoto similarity)", "Conformal interval\n(80% nominal coverage)"]
for k, s in enumerate(checks):
    yc = T3 - 3.4 - 4.4 * k
    tick(80.4, yc)
    txt(82.2, yc, s, size=8.5, box=(81.6, yc - 1.9, 98.5, yc + 1.9))

# ---------------------------------------------------------------- 4. consensus
header(43.6, 4, "Consensus: models vote and the count defines the alert")
txt(24.4, 38.8, "61 models loaded", size=9.5, weight="bold", ha="center", box=(10, 37.9, 39, 39.7))
mol_icon(6.5, 30.0, r=1.8)
txt(6.5, 24.6, "New\nmolecule", size=8.5, ha="center", box=(2, 23.2, 11, 26.0))
arrow(9.6, 30.0, 12.4, 30.0)
cols = [C_IC50] * 18 + [C_KEC] * 32 + [C_TOX] * 11
for k, c in enumerate(cols):
    rr, cc = divmod(k, 10)
    ax.add_patch(Circle((14.4 + 2.3 * cc, 36.6 - 2.3 * rr), 0.85, fc=c, ec="white", lw=0.6, zorder=4))
for x, c, lab in ((11.8, C_IC50, "18 IC50"), (19.4, C_KEC, "32 Ki, Kd, EC50"), (30.4, C_TOX, "11 Tox21")):
    ax.add_patch(Circle((x, 19.9), 0.8, fc=c, ec="none", zorder=4))
    txt(x + 1.2, 19.9, lab, size=8, box=(x, 18.9, x + 11, 20.9))
arrow(37.2, 29.0, 41.8, 29.0)
txt(54, 38.8, "Each model gives one vote", size=9.5, weight="bold", ha="center", box=(42.5, 37.9, 65.5, 39.7))
rbox(43, 33.4, 65, 37.0, NEUTRAL, ec="#9aa9c0", lw=1, r=0.9, z=3)
txt(44.0, 35.9, "Regression", size=9, weight="bold", box=(43, 33.4, 65, 37.0))
txt(44.0, 34.2, "active if predicted value ≥ 6.0", size=8.5, box=(43, 33.4, 65, 37.0))
rbox(43, 28.9, 65, 32.5, NEUTRAL, ec="#9aa9c0", lw=1, r=0.9, z=3)
txt(44.0, 31.4, "Classification", size=9, weight="bold", box=(43, 28.9, 65, 32.5))
txt(44.0, 29.7, "active if prob. ≥ model threshold", size=8.5, box=(43, 28.9, 65, 32.5))
ax.add_patch(Circle((44.2, 26.9), 0.8, fc=INK, ec="none", zorder=4))
txt(45.5, 26.9, "active vote", size=8, box=(44.8, 26.0, 53, 27.8))
ax.add_patch(Circle((54.4, 26.9), 0.8, fc="none", ec=INK, lw=1.3, zorder=4))
txt(55.7, 26.9, "inactive vote", size=8, box=(55.0, 26.0, 65.5, 27.8))
txt(43.2, 25.2, "1 vote per model: regression and\nclassification for the same endpoint\ncount as 2 votes", size=8,
    color=MUTED, va="top", ls=1.25, box=(42.8, 21.2, 65.5, 25.4))
rbox(43, 16.9, 65, 20.5, "white", ec=INK, lw=1.0, r=0.9, ls="--", z=3)
txt(44.0, 18.7, "Only deployable models whose\ndomain contains the molecule vote.", size=8, box=(43, 16.9, 65, 20.5))
arrow(65.5, 29.0, 67.6, 29.0)
txt(68, 38.8, "The count defines the alert", size=9.5, weight="bold", box=(67.5, 37.9, 98.5, 39.7))
txt(68, 37.1, "checks top to bottom; the first rule applies", size=8, color=MUTED, box=(67.5, 36.2, 98.5, 38.0))
levels = [("HIGH", HIGH, "white", "≥ 2 active votes, at least\n1 from ChEMBL and 1 from Tox21"),
          ("MODERATE", MOD, INK, "≥ 50% active votes, with\n3 or more models voting"),
          ("LOW", BAI, INK, "at least 1 active vote"),
          ("NO ALERT", NO_ALERT, INK, "no active vote")]
for i, (name, fc, tc, rule) in enumerate(levels):
    y1 = 35.4 - 4.7 * i
    y0 = y1 - 4.3
    rbox(68, y0, 79, y1, fc, r=1.0, z=3)
    txt(73.5, (y0 + y1) / 2, name, size=9, weight="bold", color=tc, ha="center", box=(68, y0, 79, y1))
    txt(80.4, (y0 + y1) / 2, rule, size=8.5, ls=1.2, box=(80.0, y0 - 0.2, 98.5, y1 + 0.2))
arrow(74.5, 17.0, 74.5, 15.9, ms=11, lw=1.8)

# result band
rbox(2, 12.2, 98, 15.6, INK, r=1.2, z=3)
txt(4, 13.9, "Result: one alert level per receptor", size=9.5, weight="bold", color="white",
    box=(3, 12.3, 41, 15.5))
for k, name in enumerate(["AR", "ERα", "ERβ", "PPARγ", "THRβ", "THRα*"]):
    chip(43.0 + 9.0 * k, 12.9, 8.0, 2.0, name, "white", tc=INK, size=8)

# footer: limits
rbox(2, 1.5, 98, 10.8, "#f1eee6", ec=RULE, lw=1.0, r=1.2, z=1)
columns = [(4.0, "Coverage",
            "Six receptors: AR, ERα, ERβ, PPARγ, THRα\nand THRβ. The insulin receptor (IR) is not yet\ncovered. *THRα has no deployable ChEMBL model."),
           (35.5, "Alert cap",
            "The high level is only possible for AR, ERα and PPARγ.\nERβ reaches moderate (no Tox21 assay) and\nTHRβ reaches low (few deployable models)."),
           (67.0, "How to read",
            "The alert sums agonist and antagonist votes and\ndoes not indicate the direction of the effect. It is a screening\nto prioritize compounds, not proof.")]
for x, cab, body_text in columns:
    txt(x, 9.2, cab, size=9.5, weight="bold", box=(x - 0.5, 8.3, x + 30.5, 10.1))
    txt(x, 7.7, body_text, size=8.5, va="top", ls=1.3, box=(x - 0.5, 2.0, x + 30.5, 7.9))


# Geometry self-check: report any text that leaves its box, leaves the figure, or overlaps another text.
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
                out.append(("outside the box", t.get_text()[:38].replace("\n", " / "),
                            tuple(round(v, 1) for v in (x0, y0, x1, y1)), box))
        if x0 < 0 or x1 > W or y0 < 0 or y1 > H:
            out.append(("outside the figure", t.get_text()[:38]))
    for i, (a, ea) in enumerate(ext):
        for b, eb in ext[i + 1:]:
            if ea[0] < eb[2] and eb[0] < ea[2] and ea[1] < eb[3] and eb[1] < ea[3]:
                out.append(("overlap", a.get_text()[:25].replace("\n", " / "), b.get_text()[:25].replace("\n", " / ")))
    return out


issues = verify()
fig.savefig("edc_model_infographic.png", dpi=220, facecolor=BG)
print(len(TEXTS), "texts |", len(issues), "issues")
for p in issues:
    print(p)
