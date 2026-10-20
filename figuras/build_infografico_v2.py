"""Infográfico da arquitetura de consenso v2 (dois eixos de evidência).

Números lidos de dados_infografico_v2.json, extraído das saídas de consenso_v2.py no painel_ref.csv.
"""
import json

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

D = json.load(open("dados_infografico_v2.json"))
MAT = np.array(D["matriz"])

W, H = 100.0, 138.0
FIGW, FIGH = 12.0, 16.56
BG, INK, MUTED, RULE = "#fbfaf7", "#1d2433", "#566070", "#d5d1c7"
FARM, TRI = "#1b3f73", "#6a4c93"            # eixo farmacológico (ChEMBL) e eixo de triagem (Tox21)
CONV, DIV, SEM, NEG = "#1b3f73", "#d9822b", "#c9cdd4", "#9fb3cd"
NEUTRO, PALE = "#e6eaf0", "#f1eee6"

fig = plt.figure(figsize=(FIGW, FIGH), facecolor=BG)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
TEXTS = []


# Draw a text item and register it for the overlap/overflow check.
def txt(x, y, s, size=8.5, weight="normal", color=INK, ha="left", va="center", style="normal", box=None, ls=1.35):
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
    txt(x0 + w / 2, y0 + h / 2, label, size=size, weight="bold", color=tc, ha="center",
        box=(x0 - 0.2, y0, x0 + w + 0.2, y0 + h))


# Numbered section header with a horizontal rule underneath.
def header(y, n, title):
    ax.add_patch(Circle((4.6, y), 2.0, fc=INK, ec="none", zorder=4))
    txt(4.6, y - 0.05, str(n), size=13, weight="bold", color="white", ha="center", box=(2.6, y - 2.0, 6.6, y + 2.0))
    txt(8.0, y, title, size=14, weight="bold", box=(7.5, y - 1.7, 98.5, y + 1.7))
    ax.plot([2, 98], [y - 2.8, y - 2.8], color=RULE, lw=1.0, zorder=1)


# Small molecule glyph standing for the input compound.
def mol_icon(cx, cy, r=1.8):
    ax.add_patch(RegularPolygon((cx, cy), 6, radius=r, orientation=0, fc="none", ec=INK, lw=1.7, zorder=5))
    ax.add_patch(Circle((cx, cy), r * 0.55, fc="none", ec=INK, lw=1.0, zorder=5))
    for ang, d, col in ((30, (1.7, 1.0), "#c0392b"), (210, (-1.7, -1.0), "#2c6fbb")):
        a = np.deg2rad(ang)
        v = (cx + r * np.cos(a), cy + r * np.sin(a))
        e = (v[0] + d[0], v[1] + d[1])
        ax.plot([v[0], e[0]], [v[1], e[1]], color=INK, lw=1.5, zorder=5, solid_capstyle="round")
        ax.add_patch(Circle(e, 0.7, fc=col, ec="white", lw=0.5, zorder=6))


# ------------------------------------------------------------------ título
txt(50, 134.3, "Consenso v2: duas linhas de evidência, lado a lado", size=20, weight="bold", ha="center",
    box=(1, 132.4, 99, 136.2))
txt(50, 130.4, "Os mesmos 61 modelos e os mesmos votos. Só a forma de resumir muda: em vez de um nível único,\n"
               "o resultado mostra o que cada tipo de dado disse e quanto dele cobriu a molécula.",
    size=10.5, color=MUTED, ha="center", box=(1, 127.8, 99, 132.4))

# ------------------------------------------------------------------ 1. o problema
header(123.4, 1, "O problema do nível único")
rbox(2, 101.5, 48.5, 119.4, "white", ec=RULE, lw=1.3, r=1.2, z=2)
txt(4.2, 117.3, "Antes: todos os votos numa só conta", size=10, weight="bold", box=(3.6, 116.1, 48, 118.5))
for k, (s, col) in enumerate((("IC50  regressão", FARM), ("IC50  classificação", FARM),
                              ("EC50  classificação", "#2b7bb9"), ("Tox21  agonista", TRI),
                              ("Tox21  antagonista", TRI))):
    y = 114.6 - 1.9 * k
    ax.add_patch(Circle((5.6, y), 0.8, fc=col, ec="none", zorder=4))
    txt(7.2, y, s, size=8.5, box=(6.6, y - 1.0, 29, y + 1.0))
txt(31.5, 110.8, "+", size=15, weight="bold", color=MUTED, ha="center", box=(30.3, 109.3, 32.7, 112.3))
chip(34.0, 109.4, 12.5, 2.8, "MODERADO", "#8a8a8a", size=8.5)
txt(4.2, 105.4, "A mesma medida conta duas vezes, os dois lados não dizem\n"
                "a mesma coisa, e quem fica fora do domínio desaparece.",
    size=8.5, va="top", ls=1.45, box=(3.6, 102.0, 48, 105.8))

rbox(51.5, 101.5, 98, 119.4, "white", ec=FARM, lw=1.5, r=1.2, z=2)
txt(53.7, 117.3, "Agora: uma linha por tipo de evidência", size=10, weight="bold", box=(53.1, 116.1, 97.5, 118.5))
for k, (nome, col, estado, cob) in enumerate(((" Farmacológico (ChEMBL)", FARM, "não avaliado", "0 de 3 endpoints no domínio"),
                                              (" Triagem (Tox21)", TRI, "positivo", "2 de 4 ensaios, antagonista"))):
    y = 113.9 - 5.6 * k
    rbox(53.7, y - 2.2, 95.8, y + 2.2, NEUTRO, ec=col, lw=1.2, r=0.9, z=3)
    ax.add_patch(Rectangle((53.7, y - 2.2), 0.9, 4.4, fc=col, ec="none", zorder=4))
    txt(55.4, y + 1.1, nome.strip(), size=9, weight="bold", color=col, box=(55.0, y + 0.1, 95.5, y + 2.1))
    txt(55.4, y - 1.1, f"{estado}  ·  {cob}", size=8.5, box=(55.0, y - 2.1, 95.5, y - 0.1))
txt(53.7, 105.4, "A triagem viu atividade; o eixo farmacológico não tinha\n"
                 "modelo capaz de opinar. O “moderado” não dizia isso.",
    size=8.5, va="top", ls=1.45, box=(53.1, 102.0, 97.5, 105.8))

# ------------------------------------------------------------------ 2. do modelo ao endpoint
header(98.4, 2, "A unidade passa a ser o endpoint, não o modelo")
mol_icon(7.0, 88.0)
txt(7.0, 83.2, "Molécula", size=8.5, ha="center", box=(2.5, 82.2, 11.5, 84.2))
arrow(10.0, 88.0, 13.0, 88.0)
txt(14.0, 94.3, "Dois modelos do mesmo endpoint", size=9.5, weight="bold", box=(13.5, 93.2, 46, 95.4))
for k, (t, sub, col) in enumerate((("Classificação", "probabilidade de ser ativo", FARM),
                                   ("Regressão", "valor de pAct e intervalo", "#2b7bb9"))):
    y = 90.6 - 5.2 * k
    rbox(14.0, y - 2.1, 45.0, y + 2.1, NEUTRO, ec="#9aa9c0", lw=1.0, r=0.9, z=3)
    txt(15.2, y + 1.0, t, size=9, weight="bold", color=col, box=(14.6, y + 0.1, 44.6, y + 2.0))
    txt(15.2, y - 1.0, sub, size=8.5, box=(14.6, y - 2.0, 44.6, y - 0.1))
txt(14.0, 81.6, "No formato antigo, estes dois viravam 2 votos.", size=8, color=MUTED, box=(13.5, 80.6, 46, 82.6))
arrow(45.8, 88.0, 49.4, 88.0)
rbox(50.5, 84.0, 98.0, 95.6, "white", ec=FARM, lw=1.5, r=1.2, z=2)
txt(52.3, 93.6, "Uma evidência por endpoint", size=10, weight="bold", color=FARM, box=(51.7, 92.5, 97.5, 94.7))
txt(52.3, 91.1, "A classificação decide se o endpoint é positivo.", size=8.5, box=(51.7, 90.0, 97.5, 92.2))
txt(52.3, 88.7, "A regressão entra como potência, com o intervalo conformal.", size=8.5, box=(51.7, 87.6, 97.5, 89.8))
txt(52.3, 86.0, "Por quê: fora do treino a classificação acerta mais (ROC-AUC 0,76 a 0,94\n"
                "nos endpoints de IC50, contra R² de −0,26 a 0,66 na regressão).",
    size=8, color=MUTED, ls=1.4, box=(51.7, 84.4, 97.5, 87.6))

# ------------------------------------------------------------------ 3. estados e cobertura
header(79.4, 3, "Cada eixo tem quatro estados, e a cobertura fica à vista")
estados = [("positivo", "ao menos um endpoint\nno domínio deu ativo", CONV, "white"),
           ("negativo", "todos os endpoints\nno domínio deram inativo", NEG, INK),
           ("não avaliado", "existe modelo, mas a molécula\nestá fora do domínio dele", SEM, INK),
           ("sem modelo", "não há modelo implantável\npara esse receptor no eixo", "#e8e4da", INK)]
for k, (nome, desc, fc, tc) in enumerate(estados):
    x0 = 2 + k * 24.5
    rbox(x0, 64.0, x0 + 22.5, 75.6, "white", ec=RULE, lw=1.2, r=1.2, z=2)
    chip(x0 + 1.4, 72.2, 19.7, 2.9, nome.upper(), fc, tc=tc, size=8.5)
    txt(x0 + 1.4, 69.4, desc, size=8.5, va="top", ls=1.45, box=(x0 + 0.9, 64.6, x0 + 22.0, 70.2))
txt(2, 61.2, "“Não avaliado” nunca é somado com “negativo”, e nunca dilui a conta do outro eixo. "
             "Esse era o viés do formato antigo:", size=8.5, box=(1.6, 60.1, 98.5, 62.3))
txt(2, 58.9, "um modelo fora do domínio desaparecia, então o mesmo sinal de triagem mudava de nível "
             "só por causa da cobertura do ChEMBL.", size=8.5, box=(1.6, 57.8, 98.5, 60.0))

# ------------------------------------------------------------------ 4. matriz de leitura
header(53.4, 4, "A leitura é o cruzamento dos dois eixos")
LX0, LY0, CW, CH = 24.0, 26.6, 17.5, 6.6     # grade 3x3
rot_col = ["positivo", "negativo", "sem cobertura"]
rot_lin = ["positivo", "negativo", "sem\ncobertura"]
LAB = [["convergente\npositivo", "divergente\nfarmacológico +", "só farmacológico\npositivo"],
       ["divergente\ntriagem +", "convergente\nnegativo", "negativo parcial"],
       ["só triagem\npositivo", "negativo parcial", "sem cobertura"]]
COR = [[CONV, DIV, "#5b8fc9"], [DIV, NEG, SEM], ["#8e74b5", SEM, SEM]]
TC = [["white", "white", "white"], ["white", INK, INK], ["white", INK, INK]]
txt(LX0 + 1.5 * CW, 49.0, "Eixo de triagem  (Tox21)", size=9.5, weight="bold", color=TRI, ha="center",
    box=(LX0, 47.9, LX0 + 3 * CW, 50.1))
for j, rc in enumerate(rot_col):
    txt(LX0 + CW * (j + 0.5), 46.6, rc, size=8.5, weight="bold", ha="center",
        box=(LX0 + CW * j, 45.5, LX0 + CW * (j + 1), 47.7))
txt(LX0 - 13.0, LY0 + 1.5 * CH, "Eixo farmacológico  (ChEMBL)", size=9.5, weight="bold", color=FARM,
    ha="center", va="center", box=None).set_rotation(90)
for i in range(3):
    for j in range(3):
        x0, y1 = LX0 + CW * j, LY0 + CH * (2 - i) + CH
        rbox(x0 + 0.35, y1 - CH + 0.35, x0 + CW - 0.35, y1 - 0.35, COR[i][j], r=0.9, z=3)
        txt((x0 + CW / 2), y1 - CH / 2 + 1.0, LAB[i][j], size=8.5, weight="bold", color=TC[i][j], ha="center",
            box=(x0 + 0.2, y1 - CH + 0.3, x0 + CW - 0.2, y1 - 0.3), ls=1.2)
        n = MAT[i][j]
        txt((x0 + CW / 2), y1 - CH / 2 - 2.0, f"{n} par" + ("" if n == 1 else "es"), size=8, color=TC[i][j],
            ha="center", box=(x0 + 0.2, y1 - CH + 0.3, x0 + CW - 0.2, y1 - 0.3))
    txt(LX0 - 1.2, LY0 + CH * (2 - i) + CH / 2, rot_lin[i], size=8.5, weight="bold", ha="right", va="center",
        box=(LX0 - 11.0, LY0 + CH * (2 - i) + 0.2, LX0 - 1.0, LY0 + CH * (3 - i) - 0.2), ls=1.2)
txt(LX0 + 3 * CW + 2.0, LY0 + 2.5 * CH, "No painel de referência:", size=8.5, weight="bold",
    box=(LX0 + 3 * CW + 1.5, LY0 + 2.5 * CH - 1.1, 98.5, LY0 + 2.5 * CH + 1.1))
txt(LX0 + 3 * CW + 2.0, LY0 + 1.55 * CH,
    f"{D['total']} pares molécula-receptor.\nO formato antigo mostrava só\n{D['n_v1']}: os outros {D['n_novos']} não\n"
    f"ganhavam linha, e a falta\nde linha era lida como\n“sem alerta”.",
    size=8.5, va="center", ls=1.5, box=(LX0 + 3 * CW + 1.5, LY0 + 0.4 * CH, 98.5, LY0 + 2.5 * CH - 1.3))
txt(2, 23.0, "Onde foram os antigos níveis:", size=9, weight="bold", box=(1.6, 21.9, 30, 24.1))
antes = [("os 10 “alto”", "viram convergente positivo: os dois eixos de acordo"),
         ("os 10 “moderado”", "eram três coisas: 5 só triagem, 4 só farmacológico, 1 divergente"),
         ("os 40 “sem alerta”", "só 1 é negativo nos dois eixos; os outros 39 têm um eixo sem cobertura")]
for k, (a, b) in enumerate(antes):
    y = 20.6 - 2.4 * k
    txt(2, y, a, size=8.5, weight="bold", box=(1.6, y - 1.1, 23.5, y + 1.1))
    txt(24.0, y, b, size=8.5, box=(23.6, y - 1.1, 98.5, y + 1.1))

# ------------------------------------------------------------------ rodapé
rbox(2, 1.6, 98, 12.2, PALE, ec=RULE, lw=1.0, r=1.2, z=1)
txt(4.0, 10.4, "Cobertura por receptor (endpoints farmacológicos / ensaios Tox21)", size=9, weight="bold",
    box=(3.5, 9.3, 50, 11.5))
for k, (rec, f_, t_) in enumerate(D["cobertura"]):
    nome = {"ESR1": "ERα", "ESR2": "ERβ", "PPARG": "PPARγ", "THRB": "THRβ", "TR": "TR"}.get(rec, rec)
    x = 4.0 + k * 7.6
    txt(x, 8.0, nome, size=8.5, weight="bold", box=(x - 0.5, 7.0, x + 7.0, 9.0))
    txt(x, 6.2, f"{f_} / {t_}", size=8.5, color=MUTED, box=(x - 0.5, 5.2, x + 7.0, 7.2))
txt(4.0, 4.0, "ERβ não tem ensaio Tox21 e THRα não tem modelo implantável: um resultado negativo neles\n"
              "vale menos do que em AR ou ERα. O rótulo TR do Tox21 não separa THRα de THRβ.",
    size=8, color=MUTED, va="top", ls=1.4, box=(3.5, 1.3, 50, 4.3))
txt(53.0, 10.4, "Limites desta versão", size=9, weight="bold", box=(52.5, 9.3, 98, 11.5))
txt(53.0, 8.6, "· Não melhora o acerto, e não foi validado por acerto: é mudança de relatório.\n"
               "· Os limiares seguem os de cada modelo (0,11 a 0,40 no Tox21 contra 0,5 no ChEMBL),\n"
               "  então as probabilidades dos dois eixos ainda não são comparáveis.\n"
               "· A direção do efeito vem só do Tox21; no eixo farmacológico fica indeterminada.\n"
               "· Não há nível único nem ordem de prioridade entre moléculas.",
    size=8, color=MUTED, va="top", ls=1.5, box=(52.5, 1.6, 98, 8.9))


# Geometry self-check: report any text that leaves its box, leaves the figure, or overlaps another text.
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
                out.append(("fora da caixa", t.get_text()[:38].replace("\n", " / "),
                            tuple(round(v, 1) for v in (x0, y0, x1, y1)), box))
        if x0 < 0 or x1 > W or y0 < 0 or y1 > H:
            out.append(("fora da figura", t.get_text()[:38]))
    for i, (a, ea) in enumerate(ext):
        for b, eb in ext[i + 1:]:
            if ea[0] < eb[2] and eb[0] < ea[2] and ea[1] < eb[3] and eb[1] < ea[3]:
                out.append(("sobreposição", a.get_text()[:25].replace("\n", " / "), b.get_text()[:25].replace("\n", " / ")))
    return out


problemas = verificar()
fig.savefig("infografico_consenso_v2.png", dpi=220, facecolor=BG)
print(len(TEXTS), "textos |", len(problemas), "problemas")
for p in problemas:
    print(p)
