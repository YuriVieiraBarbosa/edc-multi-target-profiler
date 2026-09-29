"""Predição conformal indutiva e limiar de decisão otimizado.

Substitui o domínio de aplicabilidade binário por percentil de Tanimoto por um
intervalo de predição com cobertura nominal declarada, verificada empiricamente
no teste externo por scaffold. Para classificação, escolhe o limiar por CV
(máximo MCC) em vez de 0,5 fixo.
"""
import os, sys, glob, json, warnings
sys.path.append(os.getcwd())
import numpy as np, pandas as pd, joblib
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_predict
from sklearn.metrics import matthews_corrcoef, balanced_accuracy_score, roc_auc_score
import features as F
import modeling as M

warnings.filterwarnings("ignore")
SEED, ALPHA = 42, 0.20          # cobertura nominal de 80%
os.makedirs("models_final", exist_ok=True)
os.makedirs("results", exist_ok=True)


def dificuldade(modelo, X):
    """Estimativa de dificuldade por composto: dispersão entre as árvores."""
    est = getattr(modelo, "estimators_", None)
    if est is None:
        return np.ones(len(X))
    try:
        preds = np.vstack([e.predict(X) for e in np.ravel(est)])
        return preds.std(axis=0) + 1e-6
    except Exception:
        return np.ones(len(X))


def rodar(path, perf):
    ep = os.path.basename(path).replace("dataset_", "").replace(".csv", "")
    df = pd.read_csv(path)
    X, bits, names, ok = F.featurize(df["std_smiles"].tolist())
    df = df[ok].reset_index(drop=True)
    y = df["pIC50_final"].to_numpy()
    ycls = (y >= M.ACTIVE_THRESHOLD).astype(int)
    tr, te = F.scaffold_split(df["std_smiles"].tolist(), test_size=0.2, seed=SEED)

    sub = perf[(perf.endpoint == ep) & (perf.tarefa == "regressao")]
    best = sub.loc[sub["Q2_cv"].idxmax()]
    est, _ = M.reg_zoo(X.shape[1])[best["modelo"]]
    est.set_params(**json.loads(best["melhores_params"]))

    # ICP: treino próprio / calibração
    pt, cal = train_test_split(tr, test_size=0.25, random_state=SEED)
    est.fit(X[pt], y[pt])
    r_cal = np.abs(y[cal] - est.predict(X[cal]))
    d_cal = dificuldade(est, X[cal])
    n = len(cal)
    q_lvl = min(1.0, np.ceil((n + 1) * (1 - ALPHA)) / n)
    q_abs = float(np.quantile(r_cal, q_lvl))
    q_norm = float(np.quantile(r_cal / d_cal, q_lvl))

    yhat = est.predict(X[te])
    d_te = dificuldade(est, X[te])
    cob_abs = float(np.mean(np.abs(y[te] - yhat) <= q_abs))
    cob_norm = float(np.mean(np.abs(y[te] - yhat) <= q_norm * d_te))
    lin = {"endpoint": ep, "modelo": best["modelo"], "alpha": ALPHA,
           "cobertura_nominal": 1 - ALPHA, "n_calibracao": n,
           "q_absoluto": round(q_abs, 3), "largura_intervalo_fixo": round(2 * q_abs, 3),
           "cobertura_empirica_fixo": round(cob_abs, 3),
           "cobertura_empirica_normalizado": round(cob_norm, 3),
           "largura_mediana_normalizado": round(float(2 * q_norm * np.median(d_te)), 3),
           "q_normalizado": round(q_norm, 4)}

    # limiar de classificação por CV
    subc = perf[(perf.endpoint == ep) & (perf.tarefa == "classificacao")]
    if len(subc):
        bc = subc.loc[subc["ROC_AUC_cv"].idxmax()]
        estc, _ = M.clf_zoo()[bc["modelo"]]
        estc.set_params(**json.loads(bc["melhores_params"]))
        skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
        p_cv = cross_val_predict(estc, X[tr], ycls[tr], cv=skf, method="predict_proba", n_jobs=1)[:, 1]
        grid = np.round(np.linspace(0.05, 0.95, 91), 3)
        mccs = [matthews_corrcoef(ycls[tr], (p_cv >= t).astype(int)) for t in grid]
        thr = float(grid[int(np.argmax(mccs))])
        estc.fit(X[tr], ycls[tr])
        p = estc.predict_proba(X[te])[:, 1]
        lin.update({"modelo_clf": bc["modelo"], "limiar_otimizado": thr,
                    "MCC_teste_0.5": round(matthews_corrcoef(ycls[te], (p >= .5).astype(int)), 3),
                    "MCC_teste_limiar": round(matthews_corrcoef(ycls[te], (p >= thr).astype(int)), 3),
                    "BalAcc_teste_0.5": round(balanced_accuracy_score(ycls[te], (p >= .5).astype(int)), 3),
                    "BalAcc_teste_limiar": round(balanced_accuracy_score(ycls[te], (p >= thr).astype(int)), 3)})
        # grava o limiar no modelo de produção
        fp = f"models_final/{ep}_classificacao.joblib"
        if os.path.exists(fp):
            m = joblib.load(fp); m["limiar"] = thr
            joblib.dump(m, fp, compress=3)
    # grava a calibração conformal no modelo de produção
    fp = f"models_final/{ep}_regressao.joblib"
    if os.path.exists(fp):
        m = joblib.load(fp)
        m.update({"conformal_alpha": ALPHA, "conformal_q_abs": q_abs,
                  "conformal_q_norm": q_norm, "conformal_usa_dificuldade": best["modelo"] in
                  ("RandomForest",)})
        joblib.dump(m, fp, compress=3)
    return lin


if __name__ == "__main__":
    perf = pd.read_csv("results/desempenho_modelos.csv")
    out = [rodar(p, perf) for p in sorted(glob.glob("data/dataset_*.csv"))]
    d = pd.DataFrame(out)
    d.to_csv("results/conformal_e_limiares.csv", index=False)
    print(d.to_string(index=False))
