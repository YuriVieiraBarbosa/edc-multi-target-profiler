"""Modelos das camadas Ki / Kd / EC50, mesmo protocolo dos modelos de IC50."""
import os, sys, glob, json, warnings
sys.path.append(os.getcwd())
import numpy as np, pandas as pd, joblib
from sklearn.model_selection import (GridSearchCV, KFold, StratifiedKFold,
                                     cross_val_predict, train_test_split)
from sklearn.metrics import (r2_score, mean_squared_error, mean_absolute_error,
                             roc_auc_score, average_precision_score,
                             balanced_accuracy_score, matthews_corrcoef)
import features as F
import modeling as M

warnings.filterwarnings("ignore")
SEED, ALPHA, ATIVO = 42, 0.20, 6.0
os.makedirs("models_layers", exist_ok=True)


def rodar(path):
    ep = os.path.basename(path).replace("dataset_", "").replace(".csv", "")
    df = pd.read_csv(path)
    X, bits, names, ok = F.featurize(df["std_smiles"].tolist())
    df = df[ok].reset_index(drop=True)
    y = df["pAct_final"].to_numpy()
    ycls = (y >= ATIVO).astype(int)
    tr, te = F.scaffold_split(df["std_smiles"].tolist(), test_size=0.2, seed=SEED)
    kf = KFold(5, shuffle=True, random_state=SEED)
    rows, best, bq2 = [], None, -np.inf
    for nome, (est, grid) in M.reg_zoo(X.shape[1]).items():
        nj = 1 if nome in ("RandomForest", "XGBoost") else -1
        gs = GridSearchCV(est, grid, cv=kf, scoring="r2", n_jobs=nj); gs.fit(X[tr], y[tr])
        ycv = cross_val_predict(gs.best_estimator_, X[tr], y[tr], cv=kf, n_jobs=nj)
        q2 = r2_score(y[tr], ycv)
        m = gs.best_estimator_.fit(X[tr], y[tr]); yh = m.predict(X[te])
        rows.append({"endpoint": ep, "tarefa": "regressao", "modelo": nome, "Q2_cv": q2,
                     "R2_teste": r2_score(y[te], yh),
                     "RMSE_teste": float(np.sqrt(mean_squared_error(y[te], yh))),
                     "MAE_teste": mean_absolute_error(y[te], yh),
                     "params": json.dumps(gs.best_params_)})
        if q2 > bq2: bq2, best = q2, (nome, gs.best_estimator_, json.dumps(gs.best_params_))
    nome, _, prm = best
    rng = np.random.default_rng(SEED)
    est0, _ = M.reg_zoo(X.shape[1])[nome]; est0.set_params(**json.loads(prm))
    scr = [r2_score(yp, cross_val_predict(est0, X[tr], yp, cv=kf, n_jobs=1))
           for yp in (rng.permutation(y[tr]) for _ in range(5 if len(tr) > 800 else 10))]
    # conformal
    pt, cal = train_test_split(tr, test_size=0.25, random_state=SEED)
    est0.fit(X[pt], y[pt])
    r_cal = np.abs(y[cal] - est0.predict(X[cal]))
    n = len(cal); q = float(np.quantile(r_cal, min(1.0, np.ceil((n+1)*(1-ALPHA))/n)))
    est0.fit(X[pt], y[pt]); cob = float(np.mean(np.abs(y[te] - est0.predict(X[te])) <= q))
    # produção em 100%
    estf, _ = M.reg_zoo(X.shape[1])[nome]; estf.set_params(**json.loads(prm)); estf.fit(X, y)
    nn = np.array([np.partition(r, -2)[-2] for r in M._tanimoto_matrix(bits)])
    ad = float(np.percentile(nn, 5))
    joblib.dump({"model": estf, "endpoint": ep, "task": "regression", "train_bits": bits,
                 "ad_threshold": ad, "modelo": nome, "conformal_q_abs": q,
                 "conformal_modo": "fixo", "cobertura_scaffold": cob,
                 "nao_implantavel": len(df) < 200 or bq2 < 0.3},
                f"models_layers/{ep}_regressao.joblib", compress=3)
    res = {"endpoint": ep, "n": len(df), "n_scaffolds": df["std_smiles"].map(F.murcko_scaffold).nunique(),
           "frac_ativos": float(ycls.mean()), "melhor_modelo": nome, "Q2_cv": round(bq2, 3),
           "y_scrambling_Q2_medio": round(float(np.mean(scr)), 3),
           "cobertura_conformal_scaffold": round(cob, 3), "AD_limiar_tanimoto": round(ad, 3),
           "nao_implantavel": bool(len(df) < 200 or bq2 < 0.3)}
    if 0.05 < ycls.mean() < 0.95 and min(np.bincount(ycls[tr])) >= 10:
        skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
        cb, cauc = None, -np.inf
        for nc, (est, grid) in M.clf_zoo().items():
            nj = 1 if nc in ("RandomForest", "XGBoost") else -1
            gs = GridSearchCV(est, grid, cv=skf, scoring="roc_auc", n_jobs=nj); gs.fit(X[tr], ycls[tr])
            pcv = cross_val_predict(gs.best_estimator_, X[tr], ycls[tr], cv=skf,
                                    method="predict_proba", n_jobs=1)[:, 1]
            auc = roc_auc_score(ycls[tr], pcv)
            m = gs.best_estimator_.fit(X[tr], ycls[tr]); p = m.predict_proba(X[te])[:, 1]
            rows.append({"endpoint": ep, "tarefa": "classificacao", "modelo": nc,
                         "ROC_AUC_cv": auc, "ROC_AUC_teste": roc_auc_score(ycls[te], p),
                         "PR_AUC_teste": average_precision_score(ycls[te], p),
                         "BalAcc_teste": balanced_accuracy_score(ycls[te], (p >= .5).astype(int)),
                         "MCC_teste": matthews_corrcoef(ycls[te], (p >= .5).astype(int)),
                         "params": json.dumps(gs.best_params_)})
            if auc > cauc: cauc, cb = auc, (nc, json.dumps(gs.best_params_))
        nc, prc = cb
        estc, _ = M.clf_zoo()[nc]; estc.set_params(**json.loads(prc)); estc.fit(X, ycls)
        joblib.dump({"model": estc, "endpoint": ep, "task": "classification",
                     "threshold_pIC50": ATIVO, "train_bits": bits, "ad_threshold": ad,
                     "limiar": 0.5, "modelo": nc,
                     "nao_implantavel": bool(len(df) < 200)},
                    f"models_layers/{ep}_classificacao.joblib", compress=3)
        res["melhor_classificador"] = nc
    return res, rows


if __name__ == "__main__":
    allr, allrows = [], []
    for p in sorted(glob.glob("data_layers/dataset_*.csv")):
        r, rows = rodar(p); allr.append(r); allrows.extend(rows)
        print(f"{r['endpoint']:14s} n={r['n']:5d} {r['melhor_modelo']:12s} Q2={r['Q2_cv']:.3f} "
              f"scr={r['y_scrambling_Q2_medio']:.3f} cob={r['cobertura_conformal_scaffold']:.2f}", flush=True)
    pd.DataFrame(allr).to_csv("results/resumo_modelos_camadas.csv", index=False)
    pd.DataFrame(allrows).to_csv("results/desempenho_modelos_camadas.csv", index=False)
    print("OK")
