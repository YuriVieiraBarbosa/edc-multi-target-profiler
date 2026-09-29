"""Reajusta o melhor modelo de cada endpoint em 100% dos dados curados.

As métricas de validação continuam vindo do protocolo com teste externo por
scaffold (modeling.py); estes modelos finais são os de produção.
"""
import os, sys, json, glob
sys.path.append(os.getcwd())
import numpy as np, pandas as pd, joblib
import features as F
import modeling as M

perf = pd.read_csv("results/desempenho_modelos.csv")
os.makedirs("models_final", exist_ok=True)
linhas = []
for path in sorted(glob.glob("data/dataset_*.csv")):
    ep = os.path.basename(path).replace("dataset_", "").replace(".csv", "")
    df = pd.read_csv(path)
    X, bits, names, ok = F.featurize(df["std_smiles"].tolist())
    df = df[ok].reset_index(drop=True)
    y = df["pIC50_final"].to_numpy()
    ycls = (y >= M.ACTIVE_THRESHOLD).astype(int)
    sub = perf[(perf.endpoint == ep) & (perf.tarefa == "regressao")]
    best = sub.loc[sub["Q2_cv"].idxmax()]
    est, _ = M.reg_zoo(X.shape[1])[best["modelo"]]
    est.set_params(**json.loads(best["melhores_params"]))
    est.fit(X, y)
    nn = np.sort(M._tanimoto_matrix(bits), axis=1)[:, -2]
    ad = float(np.percentile(nn, 5))
    joblib.dump({"model": est, "feature_names": names, "endpoint": ep, "task": "regression",
                 "train_bits": bits, "ad_threshold": ad, "n_train": len(df),
                 "modelo": best["modelo"], "Q2_cv": float(best["Q2_cv"]),
                 "R2_teste_scaffold": float(best["R2_teste"])},
                f"models_final/{ep}_regressao.joblib", compress=3)
    subc = perf[(perf.endpoint == ep) & (perf.tarefa == "classificacao")]
    cname = None
    if len(subc):
        bc = subc.loc[subc["ROC_AUC_cv"].idxmax()]
        cname = bc["modelo"]
        estc, _ = M.clf_zoo()[cname]
        estc.set_params(**json.loads(bc["melhores_params"]))
        estc.fit(X, ycls)
        joblib.dump({"model": estc, "feature_names": names, "endpoint": ep,
                     "task": "classification", "threshold_pIC50": M.ACTIVE_THRESHOLD,
                     "train_bits": bits, "ad_threshold": ad, "n_train": len(df),
                     "modelo": cname, "ROC_AUC_cv": float(bc["ROC_AUC_cv"]),
                     "ROC_AUC_teste_scaffold": float(bc["ROC_AUC_teste"])},
                    f"models_final/{ep}_classificacao.joblib", compress=3)
    linhas.append({"endpoint": ep, "n_treino_final": len(df),
                   "modelo_regressao": best["modelo"], "Q2_cv": round(float(best["Q2_cv"]), 3),
                   "R2_teste_scaffold": round(float(best["R2_teste"]), 3),
                   "modelo_classificacao": cname,
                   "AD_limiar_tanimoto": round(ad, 3)})
    print(linhas[-1], flush=True)
pd.DataFrame(linhas).to_csv("results/modelos_finais.csv", index=False)
