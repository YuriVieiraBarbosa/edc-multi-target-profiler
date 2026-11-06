"""Inductive conformal prediction and optimized decision threshold.

Replaces the binary applicability domain based on Tanimoto percentile with a
prediction interval with declared nominal coverage, empirically verified
on the external scaffold test. For classification, chooses the threshold via CV
(maximum MCC) instead of a fixed 0.5.
"""
import os, sys, glob, json, warnings
sys.path.append(os.getcwd())
import numpy as np, pandas as pd, joblib
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_predict
from sklearn.metrics import matthews_corrcoef, balanced_accuracy_score, roc_auc_score
import features as F
import modeling as M

warnings.filterwarnings("ignore")
SEED, ALPHA = 42, 0.20          # nominal coverage of 80%
os.makedirs("models_final", exist_ok=True)
os.makedirs("results", exist_ok=True)


# Per-compound difficulty estimate (spread across the trees), used to scale the normalised conformal interval.
def difficulty(model_name, X):
    """Per-compound difficulty estimate: dispersion across the trees."""
    est = getattr(model_name, "estimators_", None)
    if est is None:
        return np.ones(len(X))
    try:
        preds = np.vstack([e.predict(X) for e in np.ravel(est)])
        return preds.std(axis=0) + 1e-6
    except Exception:
        return np.ones(len(X))


# Fit the conformal calibration and the MCC-optimised decision threshold for one endpoint.
def run(path, perf):
    ep = os.path.basename(path).replace("dataset_", "").replace(".csv", "")
    df = pd.read_csv(path)
    X, bits, names, ok = F.featurize(df["std_smiles"].tolist())
    df = df[ok].reset_index(drop=True)
    y = df["pIC50_final"].to_numpy()
    ycls = (y >= M.ACTIVE_THRESHOLD).astype(int)
    tr, te = F.scaffold_split(df["std_smiles"].tolist(), test_size=0.2, seed=SEED)

    sub = perf[(perf.endpoint == ep) & (perf.task == "regression")]
    best = sub.loc[sub["Q2_cv"].idxmax()]
    est, _ = M.reg_zoo(X.shape[1])[best["model_name"]]
    est.set_params(**json.loads(best["best_params"]))

    # ICP: own training / calibration
    pt, cal = train_test_split(tr, test_size=0.25, random_state=SEED)
    est.fit(X[pt], y[pt])
    r_cal = np.abs(y[cal] - est.predict(X[cal]))
    d_cal = difficulty(est, X[cal])
    n = len(cal)
    q_lvl = min(1.0, np.ceil((n + 1) * (1 - ALPHA)) / n)
    q_abs = float(np.quantile(r_cal, q_lvl))
    q_norm = float(np.quantile(r_cal / d_cal, q_lvl))

    yhat = est.predict(X[te])
    d_te = difficulty(est, X[te])
    cov_abs = float(np.mean(np.abs(y[te] - yhat) <= q_abs))
    cov_norm = float(np.mean(np.abs(y[te] - yhat) <= q_norm * d_te))
    lin = {"endpoint": ep, "model_name": best["model_name"], "alpha": ALPHA,
           "coverage_nominal": 1 - ALPHA, "n_calibration": n,
           "q_absolute": round(q_abs, 3), "fixed_interval_width": round(2 * q_abs, 3),
           "coverage_empirical_fixed": round(cov_abs, 3),
           "coverage_empirical_normalized": round(cov_norm, 3),
           "normalized_median_width": round(float(2 * q_norm * np.median(d_te)), 3),
           "q_normalized": round(q_norm, 4)}

    # classification threshold via CV
    subc = perf[(perf.endpoint == ep) & (perf.task == "classification")]
    if len(subc):
        bc = subc.loc[subc["ROC_AUC_cv"].idxmax()]
        est_clf, _ = M.clf_zoo()[bc["model_name"]]
        est_clf.set_params(**json.loads(bc["best_params"]))
        skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
        p_cv = cross_val_predict(est_clf, X[tr], ycls[tr], cv=skf, method="predict_proba", n_jobs=1)[:, 1]
        grid = np.round(np.linspace(0.05, 0.95, 91), 3)
        mccs = [matthews_corrcoef(ycls[tr], (p_cv >= t).astype(int)) for t in grid]
        thr = float(grid[int(np.argmax(mccs))])
        est_clf.fit(X[tr], ycls[tr])
        p = est_clf.predict_proba(X[te])[:, 1]
        lin.update({"clf_model": bc["model_name"], "optimized_threshold": thr,
                    "MCC_test_0.5": round(matthews_corrcoef(ycls[te], (p >= .5).astype(int)), 3),
                    "MCC_test_threshold": round(matthews_corrcoef(ycls[te], (p >= thr).astype(int)), 3),
                    "BalAcc_test_0.5": round(balanced_accuracy_score(ycls[te], (p >= .5).astype(int)), 3),
                    "BalAcc_test_threshold": round(balanced_accuracy_score(ycls[te], (p >= thr).astype(int)), 3)})
        # saves the threshold to the production model
        fp = f"models_final/{ep}_classification.joblib"
        if os.path.exists(fp):
            # the keys of the saved .joblib artifacts keep their original names (limiar, modelo, nao_implantavel, ...) so the existing model files stay readable
            m = joblib.load(fp); m["limiar"] = thr
            joblib.dump(m, fp, compress=3)
    # saves the conformal calibration to the production model
    fp = f"models_final/{ep}_regression.joblib"
    if os.path.exists(fp):
        m = joblib.load(fp)
        m.update({"conformal_alpha": ALPHA, "conformal_q_abs": q_abs,
                  "conformal_q_norm": q_norm, "conformal_usa_dificuldade": best["model_name"] in
                  ("RandomForest",)})
        joblib.dump(m, fp, compress=3)
    return lin


if __name__ == "__main__":
    perf = pd.read_csv("results/model_performance.csv")
    out = [run(p, perf) for p in sorted(glob.glob("data/dataset_*.csv"))]
    d = pd.DataFrame(out)
    d.to_csv("results/conformal_and_thresholds.csv", index=False)
    print(d.to_string(index=False))
