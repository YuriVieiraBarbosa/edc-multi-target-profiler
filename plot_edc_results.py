"""Gráficos dos resultados do preditor de consenso (predict_edc2.py).

Le <out>_resumo.csv e <out>_detalhado.csv (gerados por predict_edc2.py) e produz:
  <out>_dashboard.png      painel com 4 graficos de alertas
  <out>_fracao_ativa.png   mapa de calor da fração de camadas ativas

Todos os valores usam somente camadas em domínio de aplicabilidade e implantáveis
(coluna 'conta' = no_dominio & ~nao_implantavel), a mesma regra do _resumo.csv.

Uso:
    python plot_edc_results.py --out demo
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle

BASE, MEDIO, PEQ = 10, 9, 8          # escala de fontes (titulo/eixos, legenda, ticks)
ORDEM_REC = ["AR", "ESR1", "ESR2", "PPARG", "THRA", "THRB", "TR"]
CATEGORIAS = ["n.d.", "sem alerta", "baixo", "moderado", "alto"]
NIVEIS_ALERTA = ["baixo", "moderado", "alto"]
# rampa sequencial de um unico matiz (segura para daltonismo); n.d. em cinza hachurado
COR = {"n.d.": "#D9D9D9", "sem alerta": "#EEF3F8", "baixo": "#A7C7E7",
       "moderado": "#4F8BC9", "alto": "#1B3A6B"}
COR_TEXTO = {"n.d.": "#555555", "sem alerta": "#333333", "baixo": "#1A1A1A",
             "moderado": "#FFFFFF", "alto": "#FFFFFF"}

plt.rcParams.update({
    "font.size": PEQ, "axes.titlesize": BASE, "axes.labelsize": BASE,
    "xtick.labelsize": PEQ, "ytick.labelsize": PEQ, "legend.fontsize": MEDIO,
    "figure.dpi": 100, "savefig.dpi": 300, "axes.spines.top": False,
    "axes.spines.right": False,
})


# Extract the bare alert level from the summary text; anything unrecognised becomes 'sem alerta'.
def _nivel(texto):
    texto = str(texto)
    for n in ("alto", "moderado", "baixo"):
        if texto.startswith(n):
            return n
    return "sem alerta"


# Load the predictor outputs and build the molecule x receptor matrices; None when there is no data.
def carregar(out):
    """Devolve dict com resumo, matriz de niveis e ordenacoes; None se nao houver dados."""
    f_res, f_det = f"{out}_resumo.csv", f"{out}_detalhado.csv"
    for f in (f_res, f_det):
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            print(f"AVISO: {f} ausente ou vazio; gráficos não gerados")
            return None
    try:
        res = pd.read_csv(f_res)
        det = pd.read_csv(f_det)
    except pd.errors.EmptyDataError:
        print("AVISO: nenhuma molécula avaliada; gráficos não gerados")
        return None
    if res.empty or det.empty:
        print("AVISO: nenhuma molécula avaliada; gráficos não gerados")
        return None
    mols = list(dict.fromkeys(det["nome"].astype(str)))
    res["nome"] = res["nome"].astype(str)
    presentes = set(res["receptor"])
    recs = [r for r in ORDEM_REC if r in presentes] + sorted(presentes - set(ORDEM_REC))
    sem_evidencia = [r for r in dict.fromkeys(det["receptor"]) if r not in presentes]
    res["nivel"] = res["alerta"].map(_nivel)
    cat = {(r.nome, r.receptor): r.nivel for r in res.itertuples()}
    matriz = [[cat.get((m, rc), "n.d.") for rc in recs] for m in mols]
    return {"res": res, "mols": mols, "recs": recs, "matriz": matriz,
            "sem_evidencia": sem_evidencia}


# Draw a matrix of coloured cells with a text label inside each cell.
def _celulas(ax, mols, recs, cor, texto, cor_texto, hachura=None):
    """Desenha uma matriz de celulas coloridas com texto em cada celula."""
    for i in range(len(mols)):
        for j in range(len(recs)):
            h = hachura(i, j) if hachura else None
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=cor(i, j), edgecolor="white",
                                   linewidth=1.5, hatch=h))
            ax.text(j + 0.5, i + 0.5, texto(i, j), ha="center", va="center",
                    fontsize=PEQ, color=cor_texto(i, j))
    ax.set_xlim(0, len(recs))
    ax.set_ylim(len(mols), 0)
    ax.set_xticks(np.arange(len(recs)) + 0.5)
    ax.set_xticklabels(recs)
    ax.set_yticks(np.arange(len(mols)) + 0.5)
    ax.set_yticklabels(mols)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)


# Shared figure legend with one patch per alert category.
def _legenda(fig, nota):
    handles = [Patch(facecolor=COR[c], edgecolor="#888888", linewidth=0.6,
                     hatch="////" if c == "n.d." else None, label=c) for c in CATEGORIAS]
    fig.legend(handles=handles, loc="lower center", ncol=len(CATEGORIAS), frameon=False,
               bbox_to_anchor=(0.5, 0.035), handlelength=1.6, columnspacing=1.8)
    fig.text(0.5, 0.008, nota, ha="center", va="bottom", fontsize=PEQ, color="#444444")


NOTA = ("n.d.: nenhuma camada em domínio de aplicabilidade. Alertas usam apenas camadas em "
        "domínio e implantáveis.")


# Alert dashboard: one cell per molecule x receptor, coloured by alert level.
def gerar_dashboard(out):
    d = carregar(out)
    if d is None:
        return None
    print("Gerando dashboard...")
    res, mols, recs, matriz = d["res"], d["mols"], d["recs"], d["matriz"]
    nm, nr = len(mols), len(recs)
    altura = max(9.0, 5.2 + 0.32 * nm)
    fig = plt.figure(figsize=(15, altura))
    gs = fig.add_gridspec(2, 2, left=0.11, right=0.985, top=0.925, bottom=0.115,
                          hspace=0.42, wspace=0.30)
    fig.suptitle(f"Alertas de desregulação endócrina: {nm} {'molécula' if nm == 1 else 'moléculas'}, "
                 f"{nr} {'receptor' if nr == 1 else 'receptores'}",
                 fontsize=BASE, fontweight="bold", x=0.11, ha="left", y=0.985)

    # (a) nivel de alerta por molecula x receptor
    ax = fig.add_subplot(gs[0, 0])
    _celulas(ax, mols, recs,
             cor=lambda i, j: COR[matriz[i][j]],
             texto=lambda i, j: "n.d." if matriz[i][j] == "n.d." else matriz[i][j],
             cor_texto=lambda i, j: COR_TEXTO[matriz[i][j]],
             hachura=lambda i, j: "////" if matriz[i][j] == "n.d." else None)
    ax.set_title("Nível de alerta por molécula e receptor", loc="left", pad=8)

    # (b) receptores com alerta por molecula, empilhado por nivel
    ax = fig.add_subplot(gs[0, 1])
    cont = pd.DataFrame([{n: sum(1 for x in linha if x == n) for n in NIVEIS_ALERTA}
                         for linha in matriz], index=mols)
    cont["total"] = cont.sum(axis=1)
    ordem = cont.assign(nome=cont.index).sort_values(["total", "nome"],
                                                     ascending=[False, True]).index
    cont = cont.loc[ordem]
    esq = np.zeros(len(cont))
    for n in NIVEIS_ALERTA:
        ax.barh(range(len(cont)), cont[n].values, left=esq, color=COR[n],
                edgecolor="#555555", linewidth=0.5, height=0.7)
        esq += cont[n].values
    for i, t in enumerate(cont["total"].values):
        ax.text(t + 0.08, i, str(int(t)), va="center", ha="left", fontsize=PEQ)
    ax.set_yticks(range(len(cont)))
    ax.set_yticklabels(cont.index)
    ax.invert_yaxis()
    ax.set_xlim(0, nr + 0.6)
    ax.set_xticks(range(nr + 1))
    ax.set_xlabel("Receptores com alerta (n)")
    ax.set_title("Receptores com alerta por molécula", loc="left", pad=8)

    # (c) distribuicao dos niveis entre todos os pares molecula-receptor
    ax = fig.add_subplot(gs[1, 0])
    todos = [x for linha in matriz for x in linha]
    n_cat = [todos.count(c) for c in CATEGORIAS]
    ax.bar(range(len(CATEGORIAS)), n_cat, color=[COR[c] for c in CATEGORIAS],
           edgecolor="#555555", linewidth=0.5, width=0.7,
           hatch=None)
    ax.patches[0].set_hatch("////")
    for i, v in enumerate(n_cat):
        ax.text(i, v + max(n_cat) * 0.02, str(v), ha="center", va="bottom", fontsize=PEQ)
    ax.set_xticks(range(len(CATEGORIAS)))
    ax.set_xticklabels(CATEGORIAS)
    ax.set_ylim(0, max(n_cat) * 1.15 + 0.5)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.set_ylabel("Pares molécula-receptor (n)")
    ax.set_title("Distribuição dos níveis de alerta", loc="left", pad=8)

    # (d) moleculas com alerta por receptor, empilhado por nivel
    ax = fig.add_subplot(gs[1, 1])
    cont_r = pd.DataFrame([{n: sum(1 for i in range(nm) if matriz[i][j] == n)
                            for n in NIVEIS_ALERTA} for j in range(nr)], index=recs)
    base = np.zeros(nr)
    for n in NIVEIS_ALERTA:
        ax.bar(range(nr), cont_r[n].values, bottom=base, color=COR[n],
               edgecolor="#555555", linewidth=0.5, width=0.7)
        base += cont_r[n].values
    for i, t in enumerate(base):
        ax.text(i, t + 0.08, str(int(t)), ha="center", va="bottom", fontsize=PEQ)
    ax.set_xticks(range(nr))
    ax.set_xticklabels(recs)
    ax.set_ylim(0, nm + 0.8)
    ax.set_yticks(range(0, nm + 1, max(1, nm // 10)))
    ax.set_ylabel("Moléculas com alerta (n)")
    ax.set_title("Moléculas com alerta por receptor", loc="left", pad=8)

    nota = NOTA
    if d["sem_evidencia"]:
        nota += f" Sem camada em domínio para nenhuma molécula: {', '.join(d['sem_evidencia'])}."
    _legenda(fig, nota)
    arq = f"{out}_dashboard.png"
    fig.savefig(arq, facecolor="white")
    plt.close(fig)
    print(f"   OK: {arq}")
    return arq


# Heat map of the active-vote fraction per molecule x receptor.
def gerar_fracao_ativa(out):
    d = carregar(out)
    if d is None:
        return None
    print("Gerando mapa de fração ativa...")
    res, mols, recs = d["res"], d["mols"], d["recs"]
    idx = {(r.nome, r.receptor): r for r in res.itertuples()}
    cmap = plt.get_cmap("Blues")

    # Summary row for cell (i, j), or None when that pair received no votes.
    def dado(i, j):
        return idx.get((mols[i], recs[j]))

    # Cell colour, scaled by the active-vote fraction.
    def cor(i, j):
        r = dado(i, j)
        return COR["n.d."] if r is None else cmap(0.08 + 0.92 * r.fracao_ativa)

    # Cell label: the active-vote fraction and the underlying vote counts.
    def texto(i, j):
        r = dado(i, j)
        if r is None:
            return "n.d."
        return f"{r.fracao_ativa:.2f}\n{int(r.camadas_ativas)}/{int(r.camadas_avaliadas)}"

    # Text colour chosen for contrast against the cell background.
    def cor_txt(i, j):
        r = dado(i, j)
        if r is None:
            return COR_TEXTO["n.d."]
        return "#FFFFFF" if r.fracao_ativa >= 0.55 else "#1A1A1A"

    nm, nr = len(mols), len(recs)
    fig, ax = plt.subplots(figsize=(max(8.5, 1.35 * nr + 3.5), max(6.0, 0.62 * nm + 2.6)))
    fig.subplots_adjust(left=0.17, right=0.86, top=0.90, bottom=0.17)
    _celulas(ax, mols, recs, cor, texto, cor_txt,
             hachura=lambda i, j: "////" if dado(i, j) is None else None)
    ax.set_title("Fração de camadas ativas por molécula e receptor", loc="left", pad=10)
    cax = fig.add_axes([0.89, 0.17, 0.022, 0.73])
    sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.ListedColormap(
        cmap(np.linspace(0.08, 1.0, 256))), norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("Camadas ativas / avaliadas", fontsize=BASE)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=PEQ)
    nota = ("Células: fração e (camadas ativas/avaliadas). n.d.: nenhuma camada em domínio.")
    if d["sem_evidencia"]:
        nota += f"\nSem camada em domínio para nenhuma molécula: {', '.join(d['sem_evidencia'])}."
    fig.text(0.17, 0.02, nota, ha="left", va="bottom", fontsize=PEQ, color="#444444")
    arq = f"{out}_fracao_ativa.png"
    fig.savefig(arq, facecolor="white")
    plt.close(fig)
    print(f"   OK: {arq}")
    return arq


# Generate every figure for one prediction run.
def gerar_todas(out):
    print("\n" + "=" * 60)
    print("GERANDO GRÁFICOS")
    print("=" * 60)
    if carregar(out) is None:
        return []
    arquivos = [gerar_dashboard(out), gerar_fracao_ativa(out)]
    arquivos = [a for a in arquivos if a]
    if arquivos:
        print("\nGráficos gerados:")
        for a in arquivos:
            print(f"   - {a}")
    return arquivos


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True,
                    help="prefixo usado em predict_edc2.py (--out)")
    gerar_todas(ap.parse_args().out)
