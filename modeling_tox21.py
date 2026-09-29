"""Classificadores Tox21 por ensaio (química ambiental).

Mesmo protocolo dos modelos ChEMBL: split externo por scaffold, seleção por
CV 5-fold. O limiar de decisão é otimizado por CV (máximo MCC) em vez de 0,5,
porque as classes são fortemente desbalanceadas.
"""
import os, sys, glob, json, warnings
sys.path.append(os.getcwd())
import numpy as np, pandas as pd, joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_predict
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             balanced_accuracy_score, matthews_corrcoef)
from xgboost import XGBClassifier
import features as F

warnings.filterwarnings("ignore")
SEED = 42
os.makedirs("models_tox21", exist_ok=True)
os.makedirs("results", exist_ok=True)


def zoo():
    sc = lambda e: Pipeline([("vt", VarianceThreshold(0.0)), ("sc", StandardScaler()), ("est", e)])
    return {
        "Logistic": (sc(LogisticRegression(max_iter=3000, random_state=SEED,
                                           class_weight="balanced")),
                     {"est__C": [0.1, 1.0]}),
        "RandomForest": (RandomForestClassifier(n_estimators=400, n_jobs=-1, random_state=SEED,
                                                class_weight="balanced_subsample"),
                         {"max_features": ["sqrt"]}),
        "XGBoost": (XGBClassifier(n_estimators=400, learning_rate=0.05, max_depth=6,
                                  subsample=0.8, colsample_bytree=0.6, n_jobs=-1,
                                  random_state=SEED, tree_method="hist", eval_metric="logloss"),
                    {"max_depth": [4, 6]}),
    }


def melhor_limiar(y, p):
    grid = np.unique(np.round(np.linspace(0.05, 0.95, 91), 3))
    mccs = [matthews_corrcoef(y, (p >= t).astype(int)) for t in grid]
    i = int(np.argmax(mccs))
    return float(grid[i]), float(mccs[i])


def rodar(path):
    ep = os.path.basename(path).replace("dataset_", "").replace(".csv", "")
    df = pd.read_csv(path)
    X, bits, names, ok = F.featurize(df["std_smiles"].tolist())
    df = df[ok].reset_index(drop=True)
    y = df["ativo_especifico"].astype(int).to_numpy()
    if min(np.bincount(y)) < 30:
        return None, []
    tr, te = F.scaffold_split(df["std_smiles"].tolist(), test_size=0.2, seed=SEED)
    skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
    rows, best, best_auc = [], None, -np.inf
    for nome, (est, grid) in zoo().items():
        gs = GridSearchCV(est, grid, cv=skf, scoring="roc_auc", n_jobs=1)
        gs.fit(X[tr], y[tr])
        p_cv = cross_val_predict(gs.best_estimator_, X[tr], y[tr], cv=skf,
                                 method="predict_proba", n_jobs=1)[:, 1]
        auc_cv = roc_auc_score(y[tr], p_cv)
        thr, mcc_cv = melhor_limiar(y[tr], p_cv)
        m = gs.best_estimator_.fit(X[tr], y[tr])
        p = m.predict_proba(X[te])[:, 1]
        rows.append({"endpoint": ep, "modelo": nome, "ROC_AUC_cv": auc_cv,
                     "MCC_cv_no_limiar": mcc_cv, "limiar": thr,
                     "ROC_AUC_teste": roc_auc_score(y[te], p),
                     "PR_AUC_teste": average_precision_score(y[te], p),
                     "BalAcc_teste_0.5": balanced_accuracy_score(y[te], (p >= 0.5).astype(int)),
                     "BalAcc_teste_limiar": balanced_accuracy_score(y[te], (p >= thr).astype(int)),
                     "MCC_teste_0.5": matthews_corrcoef(y[te], (p >= 0.5).astype(int)),
                     "MCC_teste_limiar": matthews_corrcoef(y[te], (p >= thr).astype(int)),
                     "params": json.dumps(gs.best_params_)})
        if auc_cv > best_auc:
            best_auc, best = auc_cv, (nome, m, thr)
    nome, m, thr = best
    sim_self = np.array([np.partition(r, -2)[-2] for r in _tan(bits[tr])])
    ad = float(np.percentile(sim_self, 5))
    # refit em 100% para produção
    mf = zoo()[nome][0]
    mf.set_params(**json.loads([r for r in rows if r["modelo"] == nome][0]["params"]))
    mf.fit(X, y)
    joblib.dump({"model": mf, "endpoint": ep, "task": "classification_tox21",
                 "limiar": thr, "train_bits": bits, "ad_threshold": ad,
                 "n_train": len(df), "frac_ativos": float(y.mean()), "modelo": nome},
                f"models_tox21/{ep}.joblib", compress=3)
    res = {"endpoint": ep, "n": len(df), "frac_ativos": float(y.mean()),
           "n_scaffolds": df["std_smiles"].map(F.murcko_scaffold).nunique(),
           "melhor_modelo": nome, "limiar": thr, "AD_limiar_tanimoto": ad}
    return res, rows


def _tan(bits, block=1024):
    b = bits.astype(np.float32); s = b.sum(1)
    out = np.empty((len(b), len(b)), dtype=np.float32)
    for i in range(0, len(b), block):
        inter = b[i:i+block] @ b.T
        out[i:i+block] = inter / np.maximum(s[i:i+block, None] + s[None, :] - inter, 1e-9)
    return out


if __name__ == "__main__":
    allr, allrows = [], []
    for p in sorted(glob.glob("data_tox21/dataset_*.csv")):
        r, rows = rodar(p)
        if r is None:
            print("pulado (poucos ativos):", p, flush=True); continue
        allr.append(r); allrows.extend(rows)
        b = max((x for x in rows), key=lambda x: x["ROC_AUC_cv"])
        print(f"{r['endpoint']:34s} n={r['n']:5d} ativos={r['frac_ativos']:.3f} "
              f"{r['melhor_modelo']:12s} AUCcv={b['ROC_AUC_cv']:.3f} AUCte={b['ROC_AUC_teste']:.3f} "
              f"MCCte={b['MCC_teste_limiar']:.3f} thr={r['limiar']:.2f}", flush=True)
    pd.DataFrame(allr).to_csv("results/resumo_modelos_tox21.csv", index=False)
    pd.DataFrame(allrows).to_csv("results/desempenho_modelos_tox21.csv", index=False)
    print("OK")
