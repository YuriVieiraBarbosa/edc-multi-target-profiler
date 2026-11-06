"""Training and validation QSAR per endpoint (pIC50 regression + classification).

Protocol: external split by Bemis-Murcko scaffold (80/20); model selection
by 5-fold CV on training set; refit on full training set; evaluation on external test;
y-scrambling; applicability domain by Tanimoto similarity.
"""
import os, sys, json, glob, pickle, warnings
sys.path.append(os.getcwd())
import numpy as np
import pandas as pd
import joblib
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.svm import SVR, SVC
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import GridSearchCV, KFold, StratifiedKFold, cross_val_predict
from sklearn.metrics import (r2_score, mean_squared_error, mean_absolute_error,
                             roc_auc_score, balanced_accuracy_score, matthews_corrcoef,
                             average_precision_score)
from xgboost import XGBRegressor, XGBClassifier
import features as F

warnings.filterwarnings("ignore")
SEED = 42
ACTIVE_THRESHOLD = 6.0          # pIC50 >= 6 (IC50 <= 1 uM)
os.makedirs("models", exist_ok=True)
os.makedirs("results", exist_ok=True)


# Candidate regressors. Ridge and SVR are wrapped in variance filtering + standardisation; the tree models are not.
def reg_zoo(n_feat):
    scaled = lambda est: Pipeline([("vt", VarianceThreshold(0.0)),
                                   ("sc", StandardScaler()), ("est", est)])
    return {
        "Ridge": (scaled(Ridge(random_state=SEED)),
                  {"est__alpha": [1.0, 10.0, 100.0]}),
        "SVR": (scaled(SVR(kernel="rbf")),
                {"est__C": [1.0, 10.0], "est__gamma": ["scale"], "est__epsilon": [0.1]}),
        "RandomForest": (RandomForestRegressor(n_estimators=500, n_jobs=-1,
                                               random_state=SEED, min_samples_leaf=1),
                         {"max_features": ["sqrt", 0.2]}),
        "XGBoost": (XGBRegressor(n_estimators=600, learning_rate=0.05, max_depth=6,
                                 subsample=0.8, colsample_bytree=0.6, n_jobs=-1,
                                 random_state=SEED, tree_method="hist"),
                    {"max_depth": [4, 6], "min_child_weight": [1, 5]}),
    }


# Candidate classifiers, following the same pipeline convention as reg_zoo.
def clf_zoo():
    scaled = lambda est: Pipeline([("vt", VarianceThreshold(0.0)),
                                   ("sc", StandardScaler()), ("est", est)])
    return {
        "Logistic": (scaled(LogisticRegression(max_iter=3000, random_state=SEED)),
                     {"est__C": [0.1, 1.0, 10.0]}),
        "SVC": (scaled(SVC(kernel="rbf", probability=True, random_state=SEED)),
                {"est__C": [1.0, 10.0]}),
        "RandomForest": (RandomForestClassifier(n_estimators=500, n_jobs=-1,
                                                random_state=SEED, class_weight="balanced_subsample"),
                         {"max_features": ["sqrt", 0.2]}),
        "XGBoost": (XGBClassifier(n_estimators=600, learning_rate=0.05, max_depth=6,
                                  subsample=0.8, colsample_bytree=0.6, n_jobs=-1,
                                  random_state=SEED, tree_method="hist",
                                  eval_metric="logloss"),
                    {"max_depth": [4, 6], "min_child_weight": [1, 5]}),
    }


# Full protocol for one endpoint: scaffold split, model selection by 5-fold CV, external test and y-scrambling.
def run_endpoint(path):
    ep = os.path.basename(path).replace("dataset_", "").replace(".csv", "")
    df = pd.read_csv(path)
    X, bits, names, ok = F.featurize(df["std_smiles"].tolist())
    df = df[ok].reset_index(drop=True)
    y = df["pIC50_final"].to_numpy()
    ycls = (y >= ACTIVE_THRESHOLD).astype(int)
    tr, te = F.scaffold_split(df["std_smiles"].tolist(), test_size=0.2, seed=SEED)
    res = {"endpoint": ep, "n_total": len(df), "n_train": len(tr), "n_test": len(te),
           "n_scaffolds": df["std_smiles"].map(F.murcko_scaffold).nunique(),
           "frac_actives": float(ycls.mean()), "pIC50_threshold": ACTIVE_THRESHOLD}

    # ---------------- regression
    kf = KFold(5, shuffle=True, random_state=SEED)
    rows = []
    best, best_q2 = None, -np.inf
    for mname, (est, grid) in reg_zoo(X.shape[1]).items():
        gs = GridSearchCV(est, grid, cv=kf, scoring="r2", n_jobs=1 if mname in ("RandomForest", "XGBoost") else -1)
        gs.fit(X[tr], y[tr])
        yhat_cv = cross_val_predict(gs.best_estimator_, X[tr], y[tr], cv=kf,
                                    n_jobs=1 if mname in ("RandomForest", "XGBoost") else -1)
        q2 = r2_score(y[tr], yhat_cv)
        m = gs.best_estimator_.fit(X[tr], y[tr])
        yhat = m.predict(X[te])
        row = {"endpoint": ep, "task": "regression", "model_name": mname,
               "Q2_cv": q2, "RMSE_cv": float(np.sqrt(mean_squared_error(y[tr], yhat_cv))),
               "R2_test": r2_score(y[te], yhat),
               "RMSE_test": float(np.sqrt(mean_squared_error(y[te], yhat))),
               "MAE_test": mean_absolute_error(y[te], yhat),
               "best_params": json.dumps(gs.best_params_)}
        rows.append(row)
        if q2 > best_q2:
            best_q2, best = q2, (mname, m, yhat, yhat_cv)

    mname, m, yhat, yhat_cv = best
    # y-scrambling on the best model
    rng = np.random.default_rng(SEED)
    scr = []
    n_scr = 10 if len(tr) < 800 else 5
    for _ in range(n_scr):
        yp = rng.permutation(y[tr])
        yh = cross_val_predict(m, X[tr], yp, cv=kf, n_jobs=1)
        scr.append(r2_score(yp, yh))
    res["best_regression"] = mname
    res["y_scrambling_Q2_mean"] = float(np.mean(scr))
    res["y_scrambling_Q2_max"] = float(np.max(scr))

    # applicability domain
    sim_tr = F.tanimoto_max_to_train(bits[tr], bits[tr])  # includes self=1
    nn_tr = np.array([np.sort(row)[-2] if len(row) > 1 else 1.0 for row in
                      _tanimoto_matrix(bits[tr])])
    ad_thr = float(np.percentile(nn_tr, 5))
    sim_te = F.tanimoto_max_to_train(bits[te], bits[tr])
    in_ad = sim_te >= ad_thr
    res["AD_threshold_tanimoto"] = ad_thr
    res["frac_test_in_domain"] = float(in_ad.mean())
    if in_ad.sum() > 5:
        res["RMSE_test_in_domain"] = float(np.sqrt(mean_squared_error(y[te][in_ad], yhat[in_ad])))
    if (~in_ad).sum() > 5:
        res["RMSE_test_outside_domain"] = float(np.sqrt(mean_squared_error(y[te][~in_ad], yhat[~in_ad])))

    joblib.dump({"model": m, "feature_names": names, "endpoint": ep, "task": "regression",
                 "train_bits": bits[tr], "ad_threshold": ad_thr,
                 "train_smiles": df["std_smiles"].iloc[tr].tolist()},
                f"models/{ep}_regression.joblib", compress=3)
    np.save(f"results/{ep}_pred_test.npy",
            np.vstack([y[te], yhat, sim_te]))

    # ---------------- classification
    skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
    cbest, cbest_auc = None, -np.inf
    if 0.05 < ycls.mean() < 0.95 and min(np.bincount(ycls[tr])) >= 10:
        for mname_c, (est, grid) in clf_zoo().items():
            gs = GridSearchCV(est, grid, cv=skf, scoring="roc_auc",
                              n_jobs=1 if mname_c in ("RandomForest", "XGBoost") else -1)
            gs.fit(X[tr], ycls[tr])
            p_cv = cross_val_predict(gs.best_estimator_, X[tr], ycls[tr], cv=skf,
                                     method="predict_proba", n_jobs=1)[:, 1]
            auc_cv = roc_auc_score(ycls[tr], p_cv)
            mc = gs.best_estimator_.fit(X[tr], ycls[tr])
            p = mc.predict_proba(X[te])[:, 1]
            pred = (p >= 0.5).astype(int)
            rows.append({"endpoint": ep, "task": "classification", "model_name": mname_c,
                         "ROC_AUC_cv": auc_cv, "ROC_AUC_test": roc_auc_score(ycls[te], p),
                         "PR_AUC_test": average_precision_score(ycls[te], p),
                         "BalAcc_test": balanced_accuracy_score(ycls[te], pred),
                         "MCC_test": matthews_corrcoef(ycls[te], pred),
                         "best_params": json.dumps(gs.best_params_)})
            if auc_cv > cbest_auc:
                cbest_auc, cbest = auc_cv, (mname_c, mc)
        res["best_classification"] = cbest[0]
        joblib.dump({"model": cbest[1], "feature_names": names, "endpoint": ep,
                     "task": "classification", "threshold_pIC50": ACTIVE_THRESHOLD,
                     "train_bits": bits[tr], "ad_threshold": ad_thr},
                    f"models/{ep}_classification.joblib", compress=3)
    return res, rows


# Pairwise Tanimoto similarity matrix of a binary fingerprint block.
def _tanimoto_matrix(bits):
    b = bits.astype(np.float32)
    inter = b @ b.T
    s = b.sum(1)
    union = s[:, None] + s[None, :] - inter
    return inter / np.maximum(union, 1e-9)


if __name__ == "__main__":
    paths = sorted(glob.glob("data/dataset_*.csv"))
    all_res, all_rows = [], []
    for p in paths:
        r, rows = run_endpoint(p)
        all_res.append(r); all_rows.extend(rows)
        print(f"{r['endpoint']:10s} n={r['n_total']:5d} best={r['best_regression']:12s} "
              f"Q2={max(x['Q2_cv'] for x in rows if x['endpoint']==r['endpoint'] and x['task']=='regression'):.3f} "
              f"AD={r['frac_test_in_domain']:.2f}", flush=True)
    pd.DataFrame(all_res).to_csv("results/model_summary.csv", index=False)
    pd.DataFrame(all_rows).to_csv("results/model_performance.csv", index=False)
    print("OK")
