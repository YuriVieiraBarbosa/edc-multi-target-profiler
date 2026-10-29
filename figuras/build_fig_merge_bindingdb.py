"""Figura do experimento de mesclagem da camada BindingDB (literatura).

Uma frase: acrescentar a fracao de literatura nao muda o desempenho no espaco quimico
do ChEMBL e o transforma no espaco quimico novo.

Dados: experimento_merge.csv, 8 endpoints x 3 sementes de divisao por scaffold, com o
teste identico nas duas condicoes e algoritmo fixo (RandomForest, 400 arvores).
Uso: python build_fig_merge_bindingdb.py <experimento_merge.csv>
"""
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# --- estilo: tres tamanhos mapeados a papel (base / anotacao / tick)
BASE, ANOT, TICK = 8, 7, 6
mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": BASE,
    "axes.titlesize": BASE, "axes.labelsize": BASE,
    "xtick.labelsize": TICK, "ytick.labelsize": TICK,
    "legend.fontsize": ANOT, "axes.linewidth": 0.7,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "savefig.facecolor": "white", "figure.facecolor": "white",
})
ATUAL, NOVO, INK, MUTED, RULE = "#8c8c8c", "#1b3f73", "#1d2433", "#5a6472", "#d9d6cf"

CAMINHO = sys.argv[1] if len(sys.argv) > 1 else "experimento_merge.csv"
T = pd.read_csv(CAMINHO)

METRICAS = [("roc_auc", "ROC-AUC", 0.5, "acaso"),
            ("pr_auc", "PR-AUC", None, None),
            ("r2", "R²", 0.0, "prevê a média")]
LINHAS = [("ChEMBL_teste", "Avaliado no espaço químico do ChEMBL"),
          ("BDB_teste", "Avaliado no espaço químico novo, do BindingDB")]

# mediana das 3 sementes por endpoint, para cada condicao
med = (T.groupby(["endpoint", "teste", "condicao"])[[m for m, *_ in METRICAS]]
         .median().reset_index())

fig, axes = plt.subplots(2, 3, figsize=(6.9, 4.6))
fig.subplots_adjust(left=0.088, right=0.975, top=0.790, bottom=0.135,
                    hspace=0.40, wspace=0.40)
TEXTS = []
# piso do eixo de R2: um endpoint chega a -3,97 e comprimiria o painel inteiro,
# entao ele sai como marcador fora de escala no piso, com o valor escrito ao lado
PISO_R2 = -1.25


# Register a text artist so the geometric check can walk every string.
def txt(obj):
    TEXTS.append(obj)
    return obj


for i, (teste, cabecalho) in enumerate(LINHAS):
    sub = med[med.teste == teste]
    for j, (col, rotulo, zero, zero_txt) in enumerate(METRICAS):
        ax = axes[i, j]
        d = sub.pivot_table(index="endpoint", columns="condicao", values=col).dropna()
        piso = PISO_R2 if col == "r2" and d.min().min() < PISO_R2 else None
        fora_escala = []
        # uma linha fina por endpoint
        for ep, r in d.iterrows():
            a, b = r.A, r.B
            if piso is not None and (a < piso or b < piso):
                ax.plot([0, 1], [max(a, piso), max(b, piso)], color=MUTED, lw=0.7,
                        alpha=0.55, zorder=2)
                for x, v in ((0, a), (1, b)):
                    if v < piso:
                        fora_escala.append(v)
                        # desloca cada triangulo para que dois pontos fora de escala no
                        # mesmo x nao caiam exatamente um sobre o outro
                        ax.plot([x + 0.05 * (len(fora_escala) - 1)], [piso], marker="v",
                                ms=4.2, color=MUTED, mec="white", mew=0.5, zorder=3)
                    else:
                        ax.plot([x], [v], marker="o", ms=2.6, color=MUTED, alpha=0.55,
                                mfc="white", mew=0.6, zorder=2)
                continue
            ax.plot([0, 1], [a, b], color=MUTED, lw=0.7, alpha=0.55,
                    marker="o", ms=2.6, mfc="white", mew=0.6, zorder=2)
        # mediana entre endpoints, destacada
        ax.plot([0, 1], [d.A.median(), d.B.median()], color=NOVO, lw=2.0,
                marker="o", ms=5.0, mfc=NOVO, mec="white", mew=0.8, zorder=4,
                solid_capstyle="round")
        if zero is not None:
            ax.axhline(zero, color=INK, lw=0.6, ls=(0, (3, 2)), zorder=1)
            txt(ax.text(0.97, zero, zero_txt, size=TICK, color=MUTED, va="top",
                        ha="right", transform=ax.get_yaxis_transform(), zorder=5))
        ax.set_xlim(-0.28, 1.28)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["atual", "+ BindingDB"] if i == 1 else [])
        if i == 1:
            for lab, cor in zip(ax.get_xticklabels(), (ATUAL, NOVO)):
                lab.set_color(cor)
                lab.set_fontweight("bold")
        txt(ax.set_ylabel(rotulo, labelpad=2))
        ax.margins(y=0.16)
        if piso is not None:
            ax.set_ylim(bottom=piso - 0.18)
            if fora_escala:
                faixa = f"{min(fora_escala):.1f} a {max(fora_escala):.1f}".replace(".", ",")
                txt(ax.text(0.5, 0.012, f"▾ {len(fora_escala)} fora de escala ({faixa})",
                            transform=ax.transAxes, ha="center", va="bottom",
                            size=TICK, color=MUTED)) if len(fora_escala) > 1 else \
                txt(ax.text(0.5, 0.012, f"▾ fora de escala ({faixa.split(' a ')[0]})",
                            transform=ax.transAxes, ha="center", va="bottom",
                            size=TICK, color=MUTED))
        ax.grid(axis="y", color=RULE, lw=0.45, zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        # variacao da mediana, a unica anotacao de resultado do painel
        dif = d.B.median() - d.A.median()
        txt(ax.annotate(f"{dif:+.2f}".replace(".", ","), xy=(0.5, 0.965),
                        xycoords="axes fraction", ha="center", va="top",
                        size=ANOT, color=NOVO if dif > 0 else ATUAL, fontweight="bold"))
    # cabecalho da faixa, ancorado no topo do painel ja posicionado
    txt(fig.text(0.012, axes[i, 0].get_position().y1 + 0.016, cabecalho, size=BASE,
                 fontweight="bold", color=INK, ha="left", va="bottom"))

txt(fig.text(0.012, 0.995, "Acrescentar a fração de literatura do BindingDB ao treino",
             size=BASE + 2.5, fontweight="bold", color=INK, ha="left", va="top"))
txt(fig.text(0.012, 0.950,
             "Mesmo conjunto de teste nas duas condições; algoritmo e hiperparâmetros fixos. "
             "Mediana de 8 endpoints × 3 sementes.",
             size=ANOT, color=MUTED, ha="left", va="top"))
txt(fig.text(0.975, 0.995, "maior = melhor", size=ANOT, color=MUTED,
             ha="right", va="top", style="italic"))

fig.legend(handles=[Line2D([], [], color=MUTED, lw=0.7, marker="o", ms=2.6, mfc="white",
                           label="um endpoint"),
                    Line2D([], [], color=NOVO, lw=2.0, marker="o", ms=5.0,
                           label="mediana dos 8 endpoints")],
           loc="lower center", bbox_to_anchor=(0.54, 0.002), ncol=2,
           frameon=False, handletextpad=0.5, columnspacing=1.6)


# Geometric check: no visible text box may overlap another or leave the figure.
def verificar():
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    # ticks fora dos limites do eixo continuam existindo como Text e nao aparecem no
    # render; conta-los produziria sobreposicao fantasma com o que estiver acima do painel
    dentro = set()
    for ax_ in fig.axes:
        y0, y1 = sorted(ax_.get_ylim())
        dentro |= {id(t) for t in ax_.get_yticklabels()
                   if not (y0 <= t.get_position()[1] <= y1)}
    caixas = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
              if t.get_text().strip() and t.get_visible() and id(t) not in dentro]
    probs = [(a.get_text()[:26], b.get_text()[:26])
             for i, (a, ba) in enumerate(caixas) for b, bb in caixas[i + 1:]
             if ba.overlaps(bb)]
    fora = [t.get_text()[:26] for t, bb in caixas
            if not (fig.bbox.x0 <= bb.x0 and bb.x1 <= fig.bbox.x1
                    and fig.bbox.y0 <= bb.y0 and bb.y1 <= fig.bbox.y1)]
    return probs, fora


probs, fora = verificar()
fig.savefig("fig_merge_bindingdb.png", dpi=300, bbox_inches="tight")
print(f"{len(caixas) if (caixas := None) else len(TEXTS)} textos registrados | "
      f"sobreposicoes: {len(probs)} | fora da figura: {len(fora)}")
for p in probs[:8]:
    print("  sobrepoe:", p)
for f in fora[:8]:
    print("  fora:", f)
