"""Ablação de camadas do consenso no painel de referência (painel_ref.csv).

Entrada: painel_detalhado.csv gerado por predict_edc2.py (modelos e caminhos padrão, rodado de dentro da pasta do pipeline;
NÃO passar --model-dirs absolutos: o código compara o nome da pasta como texto e rotularia a camada Tox21 como ChEMBL).
Recalcula o alerta com a mesma regra de predict_edc2.py (reproduz o _resumo.csv) usando subconjuntos de camadas e mede,
por molécula, o maior nível de alerta entre os receptores.
Uso: python ablacao_camadas_painel.py painel_detalhado.csv painel_ref.csv [dataset_ic50.csv] [pasta_pipeline]
"""
import sys, glob
import numpy as np, pandas as pd

det_csv, painel_csv = sys.argv[1], sys.argv[2]
ic50_csv = sys.argv[3] if len(sys.argv) > 3 else None
pipe = sys.argv[4] if len(sys.argv) > 4 else None
ORD = ["sem alerta", "baixo", "moderado", "alto"]
NIV = {k: i for i, k in enumerate(ORD)}

d = pd.read_csv(det_csv)
for c in ("conta", "ativo"):
    d[c] = d[c].astype(str).str.lower().eq("true")
cls = pd.read_csv(painel_csv).set_index("nome").classe
d["grupo"] = np.where(d.camada == "Tox21", "Tox21", d.camada.str.replace("ChEMBL ", ""))


# Recompute the alert level of every molecule x receptor pair from the per-model votes.
def alerta(v):
    v = v[v.conta].copy()
    v["farm"] = v.camada.str.startswith("ChEMBL")
    v["f"], v["t"] = v.farm & v.ativo, (~v.farm) & v.ativo
    r = v.groupby(["nome", "receptor"]).agg(aval=("endpoint", "nunique"), n_at=("ativo", "sum"), ef=("f", "sum"), et=("t", "sum")).reset_index()
    r["fr"] = r.n_at / r.aval
    r["alerta"] = np.select([(r.n_at >= 2) & (r.ef >= 1) & (r.et >= 1), (r.fr >= .5) & (r.aval >= 3), r.n_at >= 1], ["alto", "moderado", "baixo"], "sem alerta")
    return r


CFG = {"Todas as camadas (atual)": ["IC50", "Ki", "Kd", "EC50", "Tox21"], "Sem Tox21 (IC50 + Ki/Kd/EC50)": ["IC50", "Ki", "Kd", "EC50"],
       "Sem Ki/Kd/EC50 (IC50 + Tox21)": ["IC50", "Tox21"], "Só ChEMBL IC50": ["IC50"], "Só ChEMBL Ki/Kd/EC50": ["Ki", "Kd", "EC50"], "Só Tox21": ["Tox21"]}
mols = sorted(d.nome.unique())
neg_cls = ("controle negativo", "controle farmaco nao-EDC")
pos = [n for n in mols if cls[n] not in neg_cls]
neg = [n for n in mols if n not in pos]

maxlv, rows = {}, []
for nome, g in CFG.items():
    a = alerta(d[d.grupo.isin(g)])
    mx = a.assign(lv=a.alerta.map(NIV)).groupby("nome").lv.max().reindex(mols).fillna(0).astype(int)
    maxlv[nome] = mx
    for lim, lab in ((1, "≥ baixo"), (2, "≥ moderado"), (3, "alto")):
        rows.append(dict(config=nome, limiar=lab, positivas_detectadas=int((mx[pos] >= lim).sum()), n_positivas=len(pos),
                         negativas_sinalizadas=int((mx[neg] >= lim).sum()), n_negativas=len(neg)))
pd.DataFrame(rows).to_csv("ablacao_camadas_painel_resumo.csv", index=False)
pd.DataFrame({k: v.map(lambda i: ORD[i]) for k, v in maxlv.items()}).assign(classe=cls.reindex(mols)).to_csv("ablacao_camadas_painel_por_molecula.csv")

if ic50_csv and pipe:   # presença no treino (primeiro bloco do InChIKey = conectividade)
    ik = d.drop_duplicates("nome").set_index("nome").inchikey
    b14 = lambda s: {x[:14] for x in s.dropna()}
    kc = lambda fs: b14(pd.concat([pd.read_csv(f, usecols=["inchikey"]).inchikey for f in fs]))
    S = {"IC50": b14(pd.read_csv(ic50_csv, usecols=["inchikey"], low_memory=False).inchikey),
         "Ki/Kd/EC50": kc(glob.glob(pipe + "/datasets_camadas_por_endpoint/dataset_*.csv")),
         "Tox21": kc(glob.glob(pipe + "/datasets_tox21_por_ensaio/dataset_*.csv"))}
    tab = pd.DataFrame({k: [ik[n][:14] in S[k] for n in mols] for k in S}, index=mols)
    tab.assign(classe=cls.reindex(mols)).to_csv("painel_presenca_no_treino.csv")
    fora = [n for n in pos if not tab.loc[n, ["IC50", "Ki/Kd/EC50"]].any()]
    dentro = [n for n in pos if n not in fora]
    out = []
    for nome in ("Só ChEMBL IC50", "Só ChEMBL Ki/Kd/EC50", "Sem Tox21 (IC50 + Ki/Kd/EC50)"):
        for lim, lab in ((1, "≥ baixo"), (2, "≥ moderado")):
            out.append(dict(config=nome, limiar=lab, fora_do_treino_chembl=f"{int((maxlv[nome][fora] >= lim).sum())}/{len(fora)}",
                            dentro_do_treino_chembl=f"{int((maxlv[nome][dentro] >= lim).sum())}/{len(dentro)}"))
    pd.DataFrame(out).to_csv("ablacao_chembl_dentro_fora_do_treino.csv", index=False)
print("OK", len(mols), "moléculas;", len(pos), "positivas;", len(neg), "negativas")
