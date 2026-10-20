#!/usr/bin/env python
"""Predição multi-alvo de potencial desregulador endócrino a partir de SMILES.

Uso:
    python predict_edc.py --smiles "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1" "OCCOC(=O)..."
    python predict_edc.py --input moleculas.csv --smiles-col smiles --out predicoes.csv

Para cada endpoint treinado devolve pIC50 previsto, probabilidade de ser ativo
(pIC50 >= 6, i.e. IC50 <= 1 uM), a similaridade de Tanimoto ao vizinho mais
próximo do treino e a flag de domínio de aplicabilidade. Predições fora do
domínio de aplicabilidade NÃO devem ser usadas para priorização.
"""
import argparse, glob, os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import joblib
import curation
import features as F

ROTULOS = {"AR": "Receptor de andrógeno", "ESR1": "Receptor de estrógeno alfa",
           "ESR2": "Receptor de estrógeno beta", "PPARG": "PPAR-gama",
           "THRA": "Receptor de hormônio tireoidiano alfa",
           "THRB": "Receptor de hormônio tireoidiano beta"}
ENSAIO = {"B": "ligação", "F": "funcional"}


# Load the regression models of the earlier single-layer run.
def carregar_modelos(model_dir="models"):
    mods = {}
    for p in sorted(glob.glob(os.path.join(model_dir, "*_regressao.joblib"))):
        ep = os.path.basename(p).replace("_regressao.joblib", "")
        mods[ep] = {"reg": joblib.load(p)}
        pc = os.path.join(model_dir, f"{ep}_classificacao.joblib")
        if os.path.exists(pc):
            mods[ep]["clf"] = joblib.load(pc)
    return mods


# Predict pIC50 for a list of SMILES using those models.
def prever(smiles_list, model_dir="models"):
    mods = carregar_modelos(model_dir)
    if not mods:
        raise SystemExit(f"nenhum modelo encontrado em {model_dir}/")
    std = [curation.standardize(s) for s in smiles_list]
    base = pd.DataFrame({"smiles_entrada": smiles_list,
                         "std_smiles": [d["std_smiles"] for d in std],
                         "inchikey": [d["inchikey"] for d in std],
                         "motivo_rejeicao": [d["reject_reason"] for d in std]})
    val = base["std_smiles"].notna().to_numpy()
    if val.sum() == 0:
        return base
    X, bits, _, ok = F.featurize(base.loc[val, "std_smiles"].tolist())
    idx = base.index[val][ok]
    for ep, m in sorted(mods.items()):
        alvo, tipo = ep.rsplit("_", 1)
        pref = f"{alvo}_{tipo}"
        base.loc[idx, f"{pref}_pIC50_prev"] = m["reg"]["model"].predict(X).round(2)
        sim = F.tanimoto_max_to_train(bits, m["reg"]["train_bits"])
        base.loc[idx, f"{pref}_tanimoto_viz"] = sim.round(3)
        base.loc[idx, f"{pref}_no_dominio"] = sim >= m["reg"]["ad_threshold"]
        if "clf" in m:
            base.loc[idx, f"{pref}_prob_ativo"] = m["clf"]["model"].predict_proba(X)[:, 1].round(3)
    probs = [c for c in base.columns if c.endswith("_prob_ativo")]
    dom = [c.replace("_prob_ativo", "_no_dominio") for c in probs]
    ativo_dom = np.zeros(len(base))
    for pc, dc in zip(probs, dom):
        ativo_dom += ((base[pc] >= 0.5) & (base[dc] == True)).astype(int)
    base["n_endpoints_ativos_no_dominio"] = ativo_dom.astype(int)
    base["perfil_multialvo"] = np.where(base["std_smiles"].isna(), "estrutura rejeitada",
                                        np.where(ativo_dom >= 3, "multi-alvo (>=3 endpoints)",
                                        np.where(ativo_dom >= 1, "seletivo (1-2 endpoints)",
                                                 "sem alerta")))
    return base


# Command-line entry point of the first predictor, kept for reference only.
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smiles", nargs="+", help="um ou mais SMILES")
    ap.add_argument("--input", help="CSV de entrada")
    ap.add_argument("--smiles-col", default="smiles")
    ap.add_argument("--name-col", default=None)
    ap.add_argument("--model-dir", default="models")
    ap.add_argument("--out", default="predicoes_edc.csv")
    a = ap.parse_args()
    if a.input:
        df = pd.read_csv(a.input)
        smi = df[a.smiles_col].astype(str).tolist()
        nomes = df[a.name_col].tolist() if a.name_col else None
    elif a.smiles:
        smi, nomes = a.smiles, None
    else:
        ap.error("informe --smiles ou --input")
    out = prever(smi, a.model_dir)
    if nomes is not None:
        out.insert(0, "nome", nomes)
    out.to_csv(a.out, index=False)
    print(f"{len(out)} moléculas -> {a.out}")
    cols = ["nome"] if nomes is not None else ["smiles_entrada"]
    cols += ["n_endpoints_ativos_no_dominio", "perfil_multialvo"]
    print(out[[c for c in cols if c in out]].to_string(index=False))


if __name__ == "__main__":
    main()
