"""Infográfico da regra de níveis por unidade de evidência (predict_edc2.py).

Números medidos no painel_ref.csv com os modelos implantados em 27/10/2026:
ERalfa tem 12 modelos implantados que colapsam em 6 unidades; a migração dos 83 pares
molécula x receptor vem da comparação entre a regra antiga e a nova.
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
FARM, TRI = "#1b3f73", "#6a4c93"          # eixo farmacológico (ChEMBL) e de triagem (Tox21)
ALTO, MOD, BAIXO, SEM, NAV = "#1B3A6B", "#4F8BC9", "#A7C7E7", "#EEF3F8", "#EADCF0"
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


# ------------------------------------------------------------------ título
txt(50, 145.6, "Como o nível de alerta é decidido", size=20, weight="bold", ha="center",
    box=(1, 143.6, 99, 147.6))
txt(50, 141.8, "Um modelo não é uma evidência. Os votos são agrupados antes de serem contados.",
    size=11, color=MUTED, ha="center", box=(1, 140.4, 99, 143.2))

# ------------------------------------------------------------------ 1. colapso
header(136.5, 1, "Primeiro os votos repetidos viram uma evidência só")
txt(4, 131.6, "Exemplo real: receptor de estrogênio alfa (ERα), que tem 12 modelos implantados.",
    size=9, color=MUTED, box=(3.6, 130.6, 80, 132.6))

txt(20, 128.0, "12 MODELOS", size=9.5, weight="bold", ha="center", box=(9, 127.0, 31, 129.0))
txt(74, 128.0, "6 UNIDADES DE EVIDÊNCIA", size=9.5, weight="bold", ha="center",
    box=(55, 127.0, 93, 129.0))

GRUPOS = [
    (FARM, ["IC50 ligação  ·  regressão", "IC50 ligação  ·  classificação"], "ESR1_B",
     "mesmos compostos, mesmas features"),
    (FARM, ["Ki ligação  ·  regressão", "Ki ligação  ·  classificação"], "ESR1_Ki_B",
     "mesmos compostos, mesmas features"),
    (FARM, ["EC50 ligação  ·  regressão", "EC50 ligação  ·  classificação"], "ESR1_EC50_B",
     "mesmos compostos, mesmas features"),
    (FARM, ["EC50 funcional  ·  regressão", "EC50 funcional  ·  classificação"], "ESR1_EC50_F",
     "mesmos compostos, mesmas features"),
    (TRI, ["Tox21 agonista  ·  HEK293T", "Tox21 agonista  ·  VM7"], "ERα agonista",
     "mesma biblioteca, outra linhagem"),
    (TRI, ["Tox21 antagonista  ·  HEK293T", "Tox21 antagonista  ·  VM7"], "ERα antagonista",
     "mesma biblioteca, outra linhagem"),
]
Y0 = 124.0
for k, (cor, modelos, unidade, razao) in enumerate(GRUPOS):
    yc = Y0 - 4.1 * k
    for j, nome in enumerate(modelos):
        yy = yc + (0.95 if j == 0 else -0.95)
        ax.add_patch(Circle((5.6, yy), 0.72, fc=cor, ec="white", lw=0.5, zorder=5))
        txt(7.0, yy, nome, size=7.6, box=(6.4, yy - 0.95, 40, yy + 0.95))
    arrow(42.0, yc, 48.2, yc, color=cor)
    txt(44.5, yc + 1.6, "viram 1", size=6.8, color=MUTED, ha="center",
        box=(40.5, yc + 0.8, 48.5, yc + 2.4))
    rbox(49, yc - 1.75, 76, yc + 1.75, "white", ec=cor, lw=1.4, r=0.9, z=3)
    ax.add_patch(Circle((51.4, yc), 0.85, fc=cor, ec="none", zorder=5))
    txt(53.2, yc, unidade, size=8.2, weight="bold", box=(52.5, yc - 1.6, 76, yc + 1.6))
    txt(77.5, yc, razao, size=7.2, color=MUTED, box=(77.0, yc - 1.1, 98.5, yc + 1.1))

txt(4, 97.6, "Sem o agrupamento, ERα daria 12 votos — mas são 6 evidências contadas duas vezes.",
    size=9, weight="bold", box=(3.6, 96.4, 90, 98.8))

# ------------------------------------------------------------------ 2. eixos
header(92.0, 2, "Cada unidade pertence a um dos dois eixos")
for x0, cor, titulo, corpo, teto in (
    (3, FARM, "Eixo farmacológico  (ChEMBL)",
     "Afinidade e potência medidas em ensaio\n"
     "dirigido ao receptor. Uma unidade por\n"
     "endpoint: receptor × medida × tipo de ensaio.\n\n"
     "Positiva quando a CLASSIFICAÇÃO acusa.\n"
     "A regressão entra como potência prevista,\n"
     "não como um voto à parte.",
     "ERα 4 unidades  ·  ERβ 4  ·  AR 4  ·  PPARγ 5  ·  THRβ 1"),
    (51, TRI, "Eixo de triagem  (Tox21)",
     "Triagem qHTS em célula, com rótulo de\n"
     "ativo/inativo do depositante. Uma unidade\n"
     "por receptor e modo de ensaio.\n\n"
     "Positiva quando QUALQUER linhagem acusa.\n"
     "As linhagens compartilham até 99,8% da\n"
     "biblioteca, então não se corroboram.",
     "ERα 2 unidades  ·  AR 2  ·  PPARγ 2  ·  TR 1  ·  ERβ 0"),
):
    x1 = x0 + 46
    rbox(x0, 66.0, x1, 88.0, "white", ec=cor, lw=1.5, r=1.2, z=2)
    rbox(x0, 83.2, x1, 88.0, cor, r=1.2, z=3)
    txt(x0 + 1.6, 85.6, titulo, size=10, weight="bold", color="white",
        box=(x0, 83.2, x1, 88.0))
    txt(x0 + 1.6, 81.6, corpo, size=8.3, va="top", ls=1.5, box=(x0 + 0.5, 70.0, x1 - 0.5, 82.2))
    txt(x0 + 1.6, 68.6, teto, size=7.6, weight="bold", color=cor,
        box=(x0 + 0.5, 67.0, x1 - 0.5, 70.2))

# ------------------------------------------------------------------ 3. escada
header(61.0, 3, "A contagem das unidades define o nível")
NIVEIS = [
    (ALTO, "white", "ALTO", "os DOIS eixos positivos", "a evidência dirigida e a de triagem\nconcordam — é a leitura mais forte"),
    (MOD, INK, "MODERADO", "2+ unidades positivas num eixo só", "há peso, mas vindo de uma única\norigem experimental"),
    (BAIXO, INK, "BAIXO", "exatamente 1 unidade positiva", "evidência isolada, sem corroboração"),
    (SEM, INK, "SEM ALERTA", "2+ unidades avaliadas, todas negativas", "foi testado e deu negativo"),
]
for k, (fc, tc, nome, regra, nota) in enumerate(NIVEIS):
    y1 = 56.4 - 5.0 * k
    y0 = y1 - 4.5
    rbox(4, y0, 24, y1, fc, ec=RULE if fc == SEM else "none", lw=0.8, r=1.0, z=3)
    txt(14, (y0 + y1) / 2, nome, size=10.5, weight="bold", color=tc, ha="center",
        box=(4, y0, 24, y1))
    txt(25.8, (y0 + y1) / 2 + 1.25, regra, size=9, weight="bold",
        box=(25.4, (y0 + y1) / 2 + 0.2, 62, (y0 + y1) / 2 + 2.3))
    txt(25.8, (y0 + y1) / 2 - 1.3, nota, size=7.8, color=MUTED, ls=1.3,
        box=(25.4, y0 - 0.5, 62, (y0 + y1) / 2 + 0.1))

rbox(64, 34.6, 96.5, 56.4, NAV, ec="#b79fc6", lw=1.3, r=1.2, ls="--", z=3)
txt(80.2, 53.6, "NÃO AVALIADO", size=10.5, weight="bold", color="#4A2B5A", ha="center",
    box=(65, 52.4, 95.5, 54.8))
txt(80.2, 50.8, "fora da escada", size=8, color="#4A2B5A", ha="center",
    box=(65, 49.9, 95.5, 51.7))
txt(65.5, 48.2, "Nenhuma unidade foi avaliada, ou só uma\n"
                "foi — e deu negativo.\n\n"
                "Não é um grau de alerta: é a ausência de\n"
                "evidência. A regra antiga chamava isso\n"
                "de “sem alerta”, que é o oposto.\n\n"
                "Fica sem avaliação a unidade cuja molécula\n"
                "cai fora do domínio do modelo.",
    size=8.0, va="top", ls=1.42, color="#3b2448", box=(65, 35.2, 96.5, 48.8))

# ------------------------------------------------------------------ 4. migração
header(32.0, 4, "O que mudou nos 83 pares do painel de referência")
MIG = [("alto", 10, [("alto", 10)]), ("moderado", 10, [("moderado", 6), ("baixo", 4)]),
       ("baixo", 23, [("baixo", 23)]), ("sem alerta", 40, [("sem alerta", 26), ("nao avaliado", 14)])]
CORES = {"alto": ALTO, "moderado": MOD, "baixo": BAIXO, "sem alerta": SEM, "nao avaliado": NAV}
TCOR = {"alto": "white", "moderado": INK, "baixo": INK, "sem alerta": INK, "nao avaliado": "#4A2B5A"}
txt(16, 27.3, "REGRA ANTIGA", size=9, weight="bold", ha="center", box=(5, 26.4, 27, 28.2))
txt(52, 27.3, "REGRA NOVA", size=9, weight="bold", ha="center", box=(41, 26.4, 63, 28.2))
yy = 24.4
for nome, n, destinos in MIG:
    alt = max(0.6 * n / 10 + 1.6, 2.3 * len(destinos))
    rbox(5, yy - alt, 27, yy, CORES[nome], ec=RULE if nome == "sem alerta" else "none", lw=0.8,
         r=0.8, z=3)
    txt(16, yy - alt / 2, f"{nome}  ·  {n}", size=8.6, weight="bold", color=TCOR[nome], ha="center",
        box=(5, yy - alt, 27, yy))
    sub = yy
    for dn, dv in destinos:
        dalt = alt * dv / n
        arrow(27.8, sub - dalt / 2, 40.5, sub - dalt / 2, color=MUTED, lw=1.0, ms=6)
        rbox(41, sub - dalt + 0.12, 63, sub - 0.12, CORES[dn],
             ec=RULE if dn == "sem alerta" else ("#b79fc6" if dn == "nao avaliado" else "none"),
             lw=0.8, r=0.7, z=3)
        txt(52, sub - dalt / 2, f"{dn}  ·  {dv}", size=8.2, weight="bold", color=TCOR[dn],
            ha="center", box=(41, sub - dalt - 0.3, 63, sub + 0.3))
        if dn == "nao avaliado":
            txt(64.5, sub - dalt / 2, "← estes pares pareciam negativos,\nmas ninguém pôde avaliá-los",
                size=8, weight="bold", color="#4A2B5A", ls=1.3,
                box=(64.0, sub - dalt - 1.0, 98.5, sub + 1.0))
        sub -= dalt
    yy -= alt + 0.8

txt(4, 5.6, "Os 10 “alto” continuaram “alto”, todos por convergência entre os dois eixos.",
    size=9, weight="bold", box=(3.6, 4.4, 72, 6.8))
txt(4, 2.6, "Leia sempre a coluna “cobertura” junto do nível: fração positiva 1,00 com cobertura 0,17 "
            "quer dizer que uma única unidade de sete opinou.",
    size=8.3, color=MUTED, box=(3.6, 1.4, 98, 3.8))


# After rendering, flag text that leaves its box, leaves the figure, or overlaps.
def verificar():
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
                out.append(("fora da caixa", t.get_text()[:40].replace("\n", " / "),
                            tuple(round(v, 1) for v in (x0, y0, x1, y1)), box))
        if x0 < 0 or x1 > W or y0 < 0 or y1 > H:
            out.append(("fora da figura", t.get_text()[:40]))
    for i, (a, ea) in enumerate(ext):
        for b, eb in ext[i + 1:]:
            if ea[0] < eb[2] and eb[0] < ea[2] and ea[1] < eb[3] and eb[1] < ea[3]:
                out.append(("sobreposicao", a.get_text()[:28].replace("\n", " / "),
                            b.get_text()[:28].replace("\n", " / ")))
    return out


problemas = verificar()
fig.savefig("infografico_niveis_alerta.png", dpi=220, facecolor=BG)
print(len(TEXTS), "textos |", len(problemas), "problemas")
for p in problemas:
    print(p)
