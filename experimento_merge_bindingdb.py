#!/usr/bin/env python
"""Mede se acrescentar a fracao de LITERATURA do BindingDB melhora o desempenho.

Desenho, para cada endpoint, com o teste IDENTICO nas duas condicoes:

    split por scaffold 80/20 nos dados do ChEMBL   -> ChEMBL_treino, ChEMBL_teste
    split por scaffold 80/20 no acrescimo BindingDB -> BDB_treino,   BDB_teste

    A (atual)  treino = ChEMBL_treino
    B (+BDB)   treino = ChEMBL_treino + BDB_treino

    avaliado separadamente em ChEMBL_teste (espaco quimico de casa) e em
    BDB_teste (o espaco novo que o acrescimo deveria cobrir)

Algoritmo e hiperparametros ficam FIXOS entre A e B, para o que se mede ser o efeito
dos dados e nao a reselecao de modelo. Isso subestima B, que poderia escolher outro
algoritmo com mais dados, e a direcao do vies esta declarada de proposito.

Nao toca em models_*: nenhum modelo implantado e alterado.

Uso: python experimento_merge_bindingdb.py <ic50_curado.csv> <saida.csv> [sementes...]
"""
import glob
import os
import sys
import warnings

warnings.filterwarnings("ignore")
sys.path.insert(0, os.getcwd())
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.metrics import (average_precision_score, matthews_corrcoef,
                             mean_absolute_error, r2_score, roc_auc_score)

import features as F

CORTE = 6.0


# Split one index array 80/20 by Bemis-Murcko scaffold, largest groups into TRAIN.
# Uses the corrected rule (acyclics get their own group), since this is new training.
def dividir(smiles, seed, test_size=0.2):
    from collections import defaultdict
    grupos = defaultdict(list)
    for i, s in enumerate(smiles):
        sc = F.murcko_scaffold(s)
        grupos[sc if sc else ("__sem_scaffold__", i)].append(i)
    conj = sorted(grupos.values(), key=lambda x: (-len(x), x[0]))
    rng = np.random.default_rng(seed)
    alvo_tr = len(smiles) - int(round(test_size * len(smiles)))
    tr, te = [], []
    for g in conj:
        (tr if len(tr) + len(g) <= alvo_tr else te).extend(g)
    if len(te) < 0.4 * int(round(test_size * len(smiles))):   # poucos scaffolds: sorteia
        idx = rng.permutation(len(smiles))
        n = int(round(test_size * len(smiles)))
        te, tr = list(idx[:n]), list(idx[n:])
    return np.array(sorted(tr)), np.array(sorted(te))


# The ChEMBL side of one endpoint: the curated dataset already on disk.
def lado_chembl(ep, ic50_path):
    p = f"datasets_camadas_por_endpoint/dataset_{ep}.csv"
    if os.path.exists(p):
        d = pd.read_csv(p, usecols=["inchikey", "std_smiles", "pAct_final"])
    else:
        ic = pd.read_csv(ic50_path, low_memory=False,
                         usecols=["inchikey", "std_smiles", "endpoint", "pIC50_final"])
        d = (ic[ic.endpoint == ep].rename(columns={"pIC50_final": "pAct_final"})
             [["inchikey", "std_smiles", "pAct_final"]])
    return d.dropna(subset=["std_smiles", "pAct_final"]).drop_duplicates("inchikey")


# The BindingDB literature addition for the same target and measure, minus anything
# already in the ChEMBL side (the point is the compounds that are NOT there).
def lado_bdb(ep, ja):
    alvo = ep.split("_")[0]
    partes = ep.split("_")
    medida = partes[1] if partes[1] in ("EC50", "Ki", "Kd", "IC50") else "IC50"
    p = f"data_bindingdb/dataset_{alvo}_{medida}_BDB.csv"
    if not os.path.exists(p):
        return pd.DataFrame(columns=["inchikey", "std_smiles", "pAct_final"])
    d = pd.read_csv(p, low_memory=False, usecols=["inchikey", "std_smiles", "pAct_final"])
    d = d.dropna(subset=["std_smiles", "pAct_final"]).drop_duplicates("inchikey")
    return d[~d.inchikey.isin(ja)]


# Fit both conditions on one endpoint and score them on the two held-out test sets.
def rodar_endpoint(ep, ic50_path, seed):
    ch = lado_chembl(ep, ic50_path)
    bd = lado_bdb(ep, set(ch.inchikey))
    if len(ch) < 50 or len(bd) < 20:
        return []
    X_ch, _, nomes, ok_ch = F.featurize(ch.std_smiles.tolist())
    ch = ch[ok_ch].reset_index(drop=True)
    X_bd, _, _, ok_bd = F.featurize(bd.std_smiles.tolist())
    bd = bd[ok_bd].reset_index(drop=True)

    tr_ch, te_ch = dividir(ch.std_smiles.tolist(), seed)
    tr_bd, te_bd = dividir(bd.std_smiles.tolist(), seed)
    y_ch, y_bd = ch.pAct_final.to_numpy(float), bd.pAct_final.to_numpy(float)

    X = {"A": X_ch[tr_ch], "B": np.vstack([X_ch[tr_ch], X_bd[tr_bd]])}
    y = {"A": y_ch[tr_ch], "B": np.concatenate([y_ch[tr_ch], y_bd[tr_bd]])}
    testes = {"ChEMBL_teste": (X_ch[te_ch], y_ch[te_ch]),
              "BDB_teste": (X_bd[te_bd], y_bd[te_bd])}

    linhas = []
    for cond in ("A", "B"):
        ycls = (y[cond] >= CORTE).astype(int)
        reg = clone(RandomForestRegressor(n_estimators=400, random_state=0, n_jobs=4))
        reg.fit(X[cond], y[cond])
        clf = None
        if len(np.unique(ycls)) > 1:
            clf = clone(RandomForestClassifier(n_estimators=400, random_state=0, n_jobs=4))
            clf.fit(X[cond], ycls)
        for nome_t, (Xt, yt) in testes.items():
            if len(yt) < 10:
                continue
            yt_cls = (yt >= CORTE).astype(int)
            pr = reg.predict(Xt)
            r = dict(endpoint=ep, semente=seed, condicao=cond, teste=nome_t,
                     n_treino=len(y[cond]), n_teste=len(yt),
                     prevalencia_teste=round(float(yt_cls.mean()), 3),
                     r2=float(r2_score(yt, pr)), mae=float(mean_absolute_error(yt, pr)),
                     rmse=float(np.sqrt(((yt - pr) ** 2).mean())))
            if clf is not None and len(np.unique(yt_cls)) > 1:
                p1 = clf.predict_proba(Xt)[:, 1]
                r.update(roc_auc=float(roc_auc_score(yt_cls, p1)),
                         pr_auc=float(average_precision_score(yt_cls, p1)),
                         mcc=float(matthews_corrcoef(yt_cls, (p1 >= 0.5).astype(int))))
            linhas.append(r)
    return linhas


if __name__ == "__main__":
    ic50_path = sys.argv[1]
    saida = sys.argv[2] if len(sys.argv) > 2 else "experimento_merge.csv"
    sementes = [int(s) for s in sys.argv[3:]] or [0, 1, 2]
    eps = pd.read_csv("endpoints_experimento.csv").endpoint_chembl.tolist() \
        if os.path.exists("endpoints_experimento.csv") else []
    if not eps:
        raise SystemExit("endpoints_experimento.csv ausente")
    todas = []
    for ep in eps:
        for sd in sementes:
            linhas = rodar_endpoint(ep, ic50_path, sd)
            todas += linhas
            print(f"{ep:14s} semente {sd}: {len(linhas)} linhas", flush=True)
    T = pd.DataFrame(todas)
    T.to_csv(saida, index=False)
    print(f"\n{len(T)} linhas gravadas em {saida}")
