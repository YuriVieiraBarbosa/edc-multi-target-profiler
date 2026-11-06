#!/usr/bin/env python
"""Measures whether adding the LITERATURE fraction from BindingDB improves performance.

Design, for each endpoint, with the test IDENTICAL in both conditions:

    80/20 scaffold split on the ChEMBL data   -> ChEMBL_train, ChEMBL_test
    80/20 scaffold split on the BindingDB addition -> BDB_train,   BDB_test

    A (current)  train = ChEMBL_train
    B (+BDB)   train = ChEMBL_train + BDB_train

    evaluated separately on ChEMBL_test (home chemical space) and on
    BDB_test (the new space the addition should cover)

Algorithm and hyperparameters remain FIXED between A and B, so that what is measured is the effect
of the data and not model reselection. This underestimates B, which could choose another
algorithm with more data, and the direction of the bias is declared on purpose.

Does not touch models_*: no deployed model is altered.

Usage: python experiment_merge_bindingdb.py <ic50_curated.csv> <output.csv> [seeds...]
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

CUTOFF = 6.0


# Split one index array 80/20 by Bemis-Murcko scaffold, largest groups into TRAIN.
# Uses the corrected rule (acyclics get their own group), since this is new training.
def split(smiles, seed, test_size=0.2):
    from collections import defaultdict
    groups = defaultdict(list)
    for i, s in enumerate(smiles):
        sc = F.murcko_scaffold(s)
        groups[sc if sc else ("__no_scaffold__", i)].append(i)
    scaffold_groups = sorted(groups.values(), key=lambda x: (-len(x), x[0]))
    rng = np.random.default_rng(seed)
    n_train_target = len(smiles) - int(round(test_size * len(smiles)))
    tr, te = [], []
    for g in scaffold_groups:
        (tr if len(tr) + len(g) <= n_train_target else te).extend(g)
    if len(te) < 0.4 * int(round(test_size * len(smiles))):   # few scaffolds: draws randomly
        idx = rng.permutation(len(smiles))
        n = int(round(test_size * len(smiles)))
        te, tr = list(idx[:n]), list(idx[n:])
    return np.array(sorted(tr)), np.array(sorted(te))


# The ChEMBL side of one endpoint: the curated dataset already on disk.
def side_chembl(ep, ic50_path):
    p = f"datasets_layers_by_endpoint/dataset_{ep}.csv"
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
def bdb_side(ep, already):
    target = ep.split("_")[0]
    parts = ep.split("_")
    measure = parts[1] if parts[1] in ("EC50", "Ki", "Kd", "IC50") else "IC50"
    p = f"data_bindingdb/dataset_{target}_{measure}_BDB.csv"
    if not os.path.exists(p):
        return pd.DataFrame(columns=["inchikey", "std_smiles", "pAct_final"])
    d = pd.read_csv(p, low_memory=False, usecols=["inchikey", "std_smiles", "pAct_final"])
    d = d.dropna(subset=["std_smiles", "pAct_final"]).drop_duplicates("inchikey")
    return d[~d.inchikey.isin(already)]


# Fit both conditions on one endpoint and score them on the two held-out test sets.
def run_one_endpoint(ep, ic50_path, seed):
    ch = side_chembl(ep, ic50_path)
    bd = bdb_side(ep, set(ch.inchikey))
    if len(ch) < 50 or len(bd) < 20:
        return []
    X_ch, _, names, ok_ch = F.featurize(ch.std_smiles.tolist())
    ch = ch[ok_ch].reset_index(drop=True)
    X_bd, _, _, ok_bd = F.featurize(bd.std_smiles.tolist())
    bd = bd[ok_bd].reset_index(drop=True)

    tr_ch, te_ch = split(ch.std_smiles.tolist(), seed)
    tr_bd, te_bd = split(bd.std_smiles.tolist(), seed)
    y_ch, y_bd = ch.pAct_final.to_numpy(float), bd.pAct_final.to_numpy(float)

    X = {"A": X_ch[tr_ch], "B": np.vstack([X_ch[tr_ch], X_bd[tr_bd]])}
    y = {"A": y_ch[tr_ch], "B": np.concatenate([y_ch[tr_ch], y_bd[tr_bd]])}
    tests = {"ChEMBL_test": (X_ch[te_ch], y_ch[te_ch]),
              "BDB_test": (X_bd[te_bd], y_bd[te_bd])}

    rows = []
    for cond in ("A", "B"):
        ycls = (y[cond] >= CUTOFF).astype(int)
        reg = clone(RandomForestRegressor(n_estimators=400, random_state=0, n_jobs=4))
        reg.fit(X[cond], y[cond])
        clf = None
        if len(np.unique(ycls)) > 1:
            clf = clone(RandomForestClassifier(n_estimators=400, random_state=0, n_jobs=4))
            clf.fit(X[cond], ycls)
        for name_t, (Xt, yt) in tests.items():
            if len(yt) < 10:
                continue
            yt_cls = (yt >= CUTOFF).astype(int)
            pr = reg.predict(Xt)
            r = dict(endpoint=ep, random_seed=seed, condition=cond, test=name_t,
                     n_train=len(y[cond]), n_test=len(yt),
                     test_prevalence=round(float(yt_cls.mean()), 3),
                     r2=float(r2_score(yt, pr)), mae=float(mean_absolute_error(yt, pr)),
                     rmse=float(np.sqrt(((yt - pr) ** 2).mean())))
            if clf is not None and len(np.unique(yt_cls)) > 1:
                p1 = clf.predict_proba(Xt)[:, 1]
                r.update(roc_auc=float(roc_auc_score(yt_cls, p1)),
                         pr_auc=float(average_precision_score(yt_cls, p1)),
                         mcc=float(matthews_corrcoef(yt_cls, (p1 >= 0.5).astype(int))))
            rows.append(r)
    return rows


if __name__ == "__main__":
    ic50_path = sys.argv[1]
    output = sys.argv[2] if len(sys.argv) > 2 else "experiment_merge.csv"
    random_seeds = [int(s) for s in sys.argv[3:]] or [0, 1, 2]
    eps = pd.read_csv("endpoints_experiment.csv").endpoint_chembl.tolist() \
        if os.path.exists("endpoints_experiment.csv") else []
    if not eps:
        raise SystemExit("endpoints_experiment.csv missing")
    all_items = []
    for ep in eps:
        for sd in random_seeds:
            rows = run_one_endpoint(ep, ic50_path, sd)
            all_items += rows
            print(f"{ep:14s} seed {sd}: {len(rows)} lines", flush=True)
    T = pd.DataFrame(all_items)
    T.to_csv(output, index=False)
    print(f"\n{len(T)} rows written to {output}")
