"""Refits the best model for each endpoint on 100% of the curated data.

Validation metrics continue to come from the protocol with external test by
scaffold (modeling.py); these final models are the production ones.
"""
import os, sys, json, glob
sys.path.append(os.getcwd())
import numpy as np, pandas as pd, joblib
import features as F
import modeling as M

perf = pd.read_csv("results/model_performance.csv")
os.makedirs("models_final", exist_ok=True)
rows = []
for path in sorted(glob.glob("data/dataset_*.csv")):
    ep = os.path.basename(path).replace("dataset_", "").replace(".csv", "")
    df = pd.read_csv(path)
    X, bits, names, ok = F.featurize(df["std_smiles"].tolist())
    df = df[ok].reset_index(drop=True)
    y = df["pIC50_final"].to_numpy()
    ycls = (y >= M.ACTIVE_THRESHOLD).astype(int)
    sub = perf[(perf.endpoint == ep) & (perf.task == "regression")]
    best = sub.loc[sub["Q2_cv"].idxmax()]
    est, _ = M.reg_zoo(X.shape[1])[best["model_name"]]
    est.set_params(**json.loads(best["best_params"]))
    est.fit(X, y)
    nn = np.sort(M._tanimoto_matrix(bits), axis=1)[:, -2]
    ad = float(np.percentile(nn, 5))
    # the keys of the saved .joblib artifacts keep their original names (limiar, modelo, nao_implantavel, ...) so the existing model files stay readable
    joblib.dump({"model": est, "feature_names": names, "endpoint": ep, "task": "regression",
                 "train_bits": bits, "ad_threshold": ad, "n_train": len(df),
                 "modelo": best["model_name"], "Q2_cv": float(best["Q2_cv"]),
                 "R2_teste_scaffold": float(best["R2_test"])},
                f"models_final/{ep}_regression.joblib", compress=3)
    subc = perf[(perf.endpoint == ep) & (perf.task == "classification")]
    cname = None
    if len(subc):
        bc = subc.loc[subc["ROC_AUC_cv"].idxmax()]
        cname = bc["model_name"]
        est_clf, _ = M.clf_zoo()[cname]
        est_clf.set_params(**json.loads(bc["best_params"]))
        est_clf.fit(X, ycls)
        joblib.dump({"model": est_clf, "feature_names": names, "endpoint": ep,
                     "task": "classification", "threshold_pIC50": M.ACTIVE_THRESHOLD,
                     "train_bits": bits, "ad_threshold": ad, "n_train": len(df),
                     "modelo": cname, "ROC_AUC_cv": float(bc["ROC_AUC_cv"]),
                     "ROC_AUC_teste_scaffold": float(bc["ROC_AUC_test"])},
                    f"models_final/{ep}_classification.joblib", compress=3)
    rows.append({"endpoint": ep, "n_train_final": len(df),
                   "regression_model": best["model_name"], "Q2_cv": round(float(best["Q2_cv"]), 3),
                   "R2_test_scaffold": round(float(best["R2_test"]), 3),
                   "classification_model": cname,
                   "AD_threshold_tanimoto": round(ad, 3)})
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv("results/final_models.csv", index=False)
