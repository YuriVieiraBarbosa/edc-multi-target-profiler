"""Análise pareada do teste justo (lê out/fair_*.json e grava out/fair_longo.csv e out/fair_resumo_pareado.csv).

Critérios definidos antes de olhar os resultados:
  - unidade de análise = endpoint; por endpoint, mediana entre sorteios do teste;
  - sorteios usados: 11-14 (ChEMBL) e 11-13 (Tox21);
  - regressão: exclui o par (endpoint, sorteio) se R² < -5 em qualquer condição (falha numérica);
  - comparação pareada por Wilcoxon; IC 95% da mediana das diferenças por bootstrap entre endpoints (4.000 réplicas, semente 0).
Uso: cd work; python analise_teste_justo.py
"""
import glob, json
import numpy as np, pandas as pd
from scipy import stats

C = ["A_sem_grandes", "B_com_grandes", "C_tamanho_igual"]
R = []
for f in sorted(glob.glob("out/fair_*.json")):
    R += json.load(open(f))
D = pd.DataFrame(R)
D = D[~((D.grupo != "tox21") & (D.semente == 15))]
D.to_csv("out/fair_longo.csv", index=False)

reg = D[D.tarefa == "regressao"]
falha = reg.groupby(["grupo", "endpoint", "semente"]).r2.transform("min") < -5
print("linhas de regressão com falha numérica:", int(falha.sum()), "de", len(reg))
regok, cl = reg[~falha], D[D.tarefa == "classificacao"]
rng = np.random.default_rng(0)


# Paired comparison of one metric between training conditions: Wilcoxon plus a bootstrap CI of the median difference.
def pareado(d, metrica, grupos):
    m = d[d.grupo.isin(grupos)].groupby(["endpoint", "condicao"])[metrica].median().unstack().dropna(subset=C)
    out = []
    for a, b, nome in (("A_sem_grandes", "B_com_grandes", "B - A"), ("A_sem_grandes", "C_tamanho_igual", "C - A"),
                       ("C_tamanho_igual", "B_com_grandes", "B - C")):
        dlt = m[b] - m[a]
        n = len(dlt)
        p = stats.wilcoxon(dlt).pvalue if n >= 6 and (dlt != 0).any() else np.nan
        bs = [np.median(rng.choice(dlt.values, n, replace=True)) for _ in range(4000)]
        out.append(dict(comparacao=nome, n_endpoints=n, mediana_ref=m[a].median(), mediana_alt=m[b].median(),
                        delta_mediano=dlt.median(), ic95_inf=np.percentile(bs, 2.5), ic95_sup=np.percentile(bs, 97.5),
                        n_diff_positiva=int((dlt > 0).sum()), n_diff_negativa=int((dlt < 0).sum()), p_wilcoxon=p))
    return pd.DataFrame(out)


tabs = []
for nome, d, met, gr in [("RMSE (regressão)", regok, "rmse", ["ic50", "layers"]), ("MAE (regressão)", regok, "mae", ["ic50", "layers"]),
                         ("R² (regressão)", regok, "r2", ["ic50", "layers"]), ("Spearman (regressão)", regok, "spearman", ["ic50", "layers"]),
                         ("ROC-AUC (ChEMBL)", cl, "roc_auc", ["ic50", "layers"]), ("ROC-AUC (Tox21)", cl, "roc_auc", ["tox21"]),
                         ("PR-AUC (ChEMBL)", cl, "pr_auc", ["ic50", "layers"]), ("PR-AUC (Tox21)", cl, "pr_auc", ["tox21"])]:
    t = pareado(d, met, gr)
    t.insert(0, "metrica", nome)
    tabs.append(t)
pd.concat(tabs).to_csv("out/fair_resumo_pareado.csv", index=False)


# Fraction of endpoint x seed pairs in which each condition beats condition A.
def consist(d, metrica, menor_melhor):
    w = d.pivot_table(index=["endpoint", "semente"], columns="condicao", values=metrica).dropna()
    f = lambda a, b: (w[b] < w[a]) if menor_melhor else (w[b] > w[a])
    return dict(metrica=metrica, pares=len(w), B_melhor_que_A=round(f(C[0], C[1]).mean(), 3), C_melhor_que_A=round(f(C[0], C[2]).mean(), 3))


pd.DataFrame([dict(grupo="regressão", **consist(regok, "rmse", True)), dict(grupo="ChEMBL", **consist(cl[cl.grupo != "tox21"], "roc_auc", False)),
              dict(grupo="Tox21", **consist(cl[cl.grupo == "tox21"], "roc_auc", False))]).to_csv("out/fair_consistencia_sementes.csv", index=False)
print("OK")
