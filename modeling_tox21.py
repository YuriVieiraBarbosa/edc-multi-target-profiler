"""Tox21 classifiers per assay (environmental chemistry).

Same protocol as the ChEMBL models: external split by scaffold, selection by
5-fold CV. The decision threshold is optimized by CV (maximum MCC) instead of 0.5,
because the classes are strongly imbalanced.
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


# Candidate classifiers for the qHTS assays (classification only: inactive AC50 values are censored).
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


# Probability threshold that maximises the MCC over a 0.05-0.95 grid.
def best_threshold(y, p):
    grid = np.unique(np.round(np.linspace(0.05, 0.95, 91), 3))
    mccs = [matthews_corrcoef(y, (p >= t).astype(int)) for t in grid]
    i = int(np.argmax(mccs))
    return float(grid[i]), float(mccs[i])


# Train and evaluate the active/inactive classifier for one Tox21 assay.
def run(path):
    ep = os.path.basename(path).replace("dataset_", "").replace(".csv", "")
    df = pd.read_csv(path)
    X, bits, names, ok = F.featurize(df["std_smiles"].tolist())
    df = df[ok].reset_index(drop=True)
    y = df["specific_active"].astype(int).to_numpy()
    if min(np.bincount(y)) < 30:
        return None, []
    tr, te = F.scaffold_split(df["std_smiles"].tolist(), test_size=0.2, seed=SEED)
    skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
    rows, best, best_auc = [], None, -np.inf
    for name, (est, grid) in zoo().items():
        gs = GridSearchCV(est, grid, cv=skf, scoring="roc_auc", n_jobs=1)
        gs.fit(X[tr], y[tr])
        p_cv = cross_val_predict(gs.best_estimator_, X[tr], y[tr], cv=skf,
                                 method="predict_proba", n_jobs=1)[:, 1]
        auc_cv = roc_auc_score(y[tr], p_cv)
        thr, mcc_cv = best_threshold(y[tr], p_cv)
        m = gs.best_estimator_.fit(X[tr], y[tr])
        p = m.predict_proba(X[te])[:, 1]
        rows.append({"endpoint": ep, "model_name": name, "ROC_AUC_cv": auc_cv,
                     "MCC_cv_at_threshold": mcc_cv, "threshold": thr,
                     "ROC_AUC_test": roc_auc_score(y[te], p),
                     "PR_AUC_test": average_precision_score(y[te], p),
                     "BalAcc_test_0.5": balanced_accuracy_score(y[te], (p >= 0.5).astype(int)),
                     "BalAcc_test_threshold": balanced_accuracy_score(y[te], (p >= thr).astype(int)),
                     "MCC_test_0.5": matthews_corrcoef(y[te], (p >= 0.5).astype(int)),
                     "MCC_test_threshold": matthews_corrcoef(y[te], (p >= thr).astype(int)),
                     "params": json.dumps(gs.best_params_)})
        if auc_cv > best_auc:
            best_auc, best = auc_cv, (name, m, thr)
    name, m, thr = best
    sim_self = np.array([np.partition(r, -2)[-2] for r in _tan(bits[tr])])
    ad = float(np.percentile(sim_self, 5))
    # refit at 100% for production
    mf = zoo()[name][0]
    mf.set_params(**json.loads([r for r in rows if r["model_name"] == name][0]["params"]))
    mf.fit(X, y)
    # the keys of the saved .joblib artifacts keep their original names (limiar, modelo, nao_implantavel, ...) so the existing model files stay readable
    joblib.dump({"model": mf, "endpoint": ep, "task": "classification_tox21",
                 "limiar": thr, "train_bits": bits, "ad_threshold": ad,
                 "n_train": len(df), "frac_ativos": float(y.mean()), "modelo": name},
                f"models_tox21/{ep}.joblib", compress=3)
    res = {"endpoint": ep, "n": len(df), "frac_actives": float(y.mean()),
           "n_scaffolds": df["std_smiles"].map(F.murcko_scaffold).nunique(),
           "best_model": name, "threshold": thr, "AD_threshold_tanimoto": ad}
    return res, rows


# Pairwise Tanimoto similarity matrix, computed in blocks to bound memory.
def _tan(bits, block=1024):
    b = bits.astype(np.float32); s = b.sum(1)
    out = np.empty((len(b), len(b)), dtype=np.float32)
    for i in range(0, len(b), block):
        inter = b[i:i+block] @ b.T
        out[i:i+block] = inter / np.maximum(s[i:i+block, None] + s[None, :] - inter, 1e-9)
    return out


# Locate the folder holding the curated assay datasets. pipeline_tox21.py writes to
# data_tox21/, but the delivered tree renamed it, so accept either name instead of failing
# with an empty glob that looks like "no datasets".
def data_folder(candidate_paths, label):
    for c in candidate_paths:
        if glob.glob(os.path.join(c, "dataset_*.csv")):
            return c
    raise SystemExit(f"no {label} dataset found. Searched for dataset_*.csv in: "
                     + ", ".join(candidate_paths) + ".\nRun the corresponding pipeline_* first, or "
                     "create a symbolic link with one of these names.")


if __name__ == "__main__":
    folder = data_folder(["data_tox21", "datasets_tox21_by_assay"], "Tox21 assays")
    print(f"datasets read from: {folder}/", flush=True)
    allr, allrows = [], []
    for p in sorted(glob.glob(f"{folder}/dataset_*.csv")):
        r, rows = run(p)
        if r is None:
            print("skipped (few actives):", p, flush=True); continue
        allr.append(r); allrows.extend(rows)
        b = max((x for x in rows), key=lambda x: x["ROC_AUC_cv"])
        print(f"{r['endpoint']:34s} n={r['n']:5d} actives={r['frac_actives']:.3f} "
              f"{r['best_model']:12s} AUCcv={b['ROC_AUC_cv']:.3f} AUCte={b['ROC_AUC_test']:.3f} "
              f"MCCte={b['MCC_test_threshold']:.3f} thr={r['threshold']:.2f}", flush=True)
    pd.DataFrame(allr).to_csv("results/model_summary_tox21.csv", index=False)
    pd.DataFrame(allrows).to_csv("results/model_performance_tox21.csv", index=False)
    print("OK")
