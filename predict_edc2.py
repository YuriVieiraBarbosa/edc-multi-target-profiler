#!/usr/bin/env python
"""Preditor de consenso multi-camada para desregulação endócrina.

Camadas independentes, nunca fundidas:
  - ChEMBL IC50   (models_final/)   regressão de pIC50 + classificação
  - ChEMBL Ki/EC50 (models_layers/) regressão de pKi / pEC50 + classificação
  - Tox21 qHTS     (models_tox21/)  classificação por ensaio, química ambiental

Uso:
    python predict_edc2.py --smiles "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1"
    python predict_edc2.py --input moleculas.csv --smiles-col smiles --name-col nome
"""
import argparse, glob, os, sys, json
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import joblib
import curation
import features as F

RECEPTOR = {"AR": "AR", "ESR1": "ERalfa", "ESR2": "ERbeta", "PPARG": "PPARgama",
            "THRA": "THRalfa", "THRB": "THRbeta", "TR": "THR (nao resolve alfa/beta)"}
ENSAIO = {"B": "ligacao", "F": "funcional"}


# Load every .joblib model and tag it with its layer, receptor and task.
def carregar(model_dirs=("models_final", "models_layers", "models_tox21")):
    camadas = []
    for d in model_dirs:
        for p in sorted(glob.glob(os.path.join(d, "*.joblib"))):
            m = joblib.load(p)
            base = os.path.basename(p)[:-7]
            if d == "models_tox21":
                _, rec, modo, cel = base.split("_", 3)
                camadas.append({"camada": "Tox21", "receptor": rec, "detalhe": f"{modo}/{cel}",
                                "tarefa": "classificacao", "modelo": m, "chave": base})
            else:
                partes = base.split("_")
                tarefa = partes[-1]
                if d == "models_final":
                    rec, atype, medida = partes[0], partes[1], "IC50"
                else:
                    rec, medida, atype = partes[0], partes[1], partes[2]
                camadas.append({"camada": f"ChEMBL {medida}", "receptor": rec,
                                "detalhe": ENSAIO.get(atype, atype),
                                "tarefa": "regressao" if tarefa == "regressao" else "classificacao",
                                "modelo": m, "chave": base})
    return camadas


# Per-compound difficulty used to scale the normalised conformal interval.
def _dificuldade(est, X):
    sub = getattr(est, "estimators_", None)
    if sub is None:
        return np.ones(len(X))
    try:
        return np.vstack([e.predict(X) for e in np.ravel(sub)]).std(axis=0) + 1e-6
    except Exception:
        return np.ones(len(X))


# Score the molecules with all models and aggregate the per-model votes for each receptor.
def prever(smiles_list, nomes=None, model_dirs=("models_final", "models_layers", "models_tox21")):
    camadas = carregar(model_dirs)
    if not camadas:
        raise SystemExit("nenhum modelo encontrado")
    std = [curation.standardize(s) for s in smiles_list]
    base = pd.DataFrame({"nome": nomes if nomes is not None else smiles_list,
                         "smiles_entrada": smiles_list,
                         "std_smiles": [d["std_smiles"] for d in std],
                         "inchikey": [d["inchikey"] for d in std],
                         "motivo_rejeicao": [d["reject_reason"] for d in std]})
    val = base["std_smiles"].notna().to_numpy()
    if val.sum() == 0:
        return base, pd.DataFrame()
    X, bits, _, ok = F.featurize(base.loc[val, "std_smiles"].tolist())
    idx = base.index[val][ok]

    linhas = []
    for c in camadas:
        m = c["modelo"]
        sim = F.tanimoto_max_to_train(bits, m["train_bits"])
        dentro = sim >= m["ad_threshold"]
        d = pd.DataFrame({"idx": idx, "camada": c["camada"], "receptor": c["receptor"],
                          "detalhe": c["detalhe"], "tarefa": c["tarefa"],
                          "endpoint": c["chave"], "tanimoto_viz": sim.round(3),
                          "no_dominio": dentro,
                          "nao_implantavel": bool(m.get("nao_implantavel", False))})
        if c["tarefa"] == "regressao":
            yh = m["model"].predict(X)
            d["valor_previsto"] = yh.round(2)
            q = m.get("conformal_q_abs")
            if q is not None:
                if m.get("conformal_modo") == "normalizado":
                    q = m["conformal_q_norm"] * _dificuldade(m["model"], X)
                d["intervalo_inf"] = (yh - q).round(2)
                d["intervalo_sup"] = (yh + q).round(2)
                d["cobertura_scaffold"] = m.get("cobertura_scaffold")
            d["ativo"] = d["valor_previsto"] >= 6.0
        else:
            p = m["model"].predict_proba(X)[:, 1]
            thr = float(m.get("limiar", 0.5))
            d["prob_ativo"] = p.round(3)
            d["limiar"] = thr
            d["ativo"] = p >= thr
        linhas.append(d)
    det = pd.concat(linhas, ignore_index=True)
    det = det.merge(base[["nome", "smiles_entrada", "inchikey"]].reset_index()
                    .rename(columns={"index": "idx"}), on="idx", how="left")

    # Usable evidence: in-domain prediction from a deployable model.
    # One model = one vote. Regression and classification of the same endpoint are two
    # independent votes, since `endpoint` already carries the task suffix.
    det["conta"] = det["no_dominio"] & ~det["nao_implantavel"]
    ev = det[det["conta"]].copy()
    ev["farmacologica"] = ev["camada"].str.startswith("ChEMBL")
    resumo = (ev.groupby(["nome", "receptor"])
                .agg(camadas_avaliadas=("endpoint", "nunique"),
                     camadas_ativas=("ativo", "sum"),
                     evidencia_farmacologica=("farmacologica", lambda s: int((s & ev.loc[s.index, "ativo"]).sum())),
                     evidencia_tox21=("farmacologica", lambda s: int((~s & ev.loc[s.index, "ativo"]).sum())))
                .reset_index())
    # Deployable models per receptor: the ceiling of votes, identical for every molecule.
    # A model left out by the domain filter is NOT counted as an inactive vote; it simply
    # does not vote, which is why the two denominators below differ.
    teto = (det[~det["nao_implantavel"]].groupby("receptor").endpoint.nunique()
            .rename("modelos_implantaveis_receptor").reset_index())
    resumo = resumo.merge(teto, on="receptor", how="left")
    # Two readings of "votes received / number of models", kept side by side:
    #   fracao_ativa              -> denominator = models that actually voted (previous scheme)
    #   fracao_sobre_implantaveis -> denominator = all deployable models of the receptor
    # The gap between them is the applicability-domain coverage for this molecule.
    resumo["fracao_ativa"] = (resumo.camadas_ativas / resumo.camadas_avaliadas).round(2)
    resumo["fracao_sobre_implantaveis"] = (resumo.camadas_ativas
                                           / resumo.modelos_implantaveis_receptor).round(2)
    resumo["cobertura_dominio"] = (resumo.camadas_avaliadas
                                   / resumo.modelos_implantaveis_receptor).round(2)
    resumo["votos"] = (resumo.camadas_ativas.astype(str) + "/" + resumo.camadas_avaliadas.astype(str)
                       + " votantes (" + resumo.camadas_ativas.astype(str) + "/"
                       + resumo.modelos_implantaveis_receptor.astype(str) + " implantaveis)")
    resumo["alerta"] = np.select(
        [(resumo.camadas_ativas >= 2) & (resumo.evidencia_farmacologica >= 1) & (resumo.evidencia_tox21 >= 1),
         (resumo.fracao_ativa >= 0.5) & (resumo.camadas_avaliadas >= 3),
         resumo.camadas_ativas >= 1],
        ["alto (farmacologico + ambiental)", "moderado (maioria das camadas)", "baixo (evidencia isolada)"],
        default="sem alerta")
    return det, resumo


# Command-line entry point of the multi-layer consensus predictor.
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smiles", nargs="+")
    ap.add_argument("--input")
    ap.add_argument("--smiles-col", default="smiles")
    ap.add_argument("--name-col", default=None)
    ap.add_argument("--model-dirs", nargs="+",
                    default=["models_final", "models_layers", "models_tox21"])
    ap.add_argument("--out", default="predicoes_consenso")
    a = ap.parse_args()
    if a.input:
        df = pd.read_csv(a.input)
        smi = df[a.smiles_col].astype(str).tolist()
        nomes = df[a.name_col].tolist() if a.name_col else None
    elif a.smiles:
        smi, nomes = a.smiles, None
    else:
        ap.error("informe --smiles ou --input")
    det, res = prever(smi, nomes, tuple(a.model_dirs))
    det.to_csv(f"{a.out}_detalhado.csv", index=False)
    res.to_csv(f"{a.out}_resumo.csv", index=False)
    print(f"{len(smi)} moléculas -> {a.out}_detalhado.csv / {a.out}_resumo.csv")
    if len(res):
        piv = res.pivot(index="nome", columns="receptor", values="alerta")
        print(piv.fillna("-").to_string())


if __name__ == "__main__":
    main()
