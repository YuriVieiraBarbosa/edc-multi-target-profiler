"""Eixo de triagem do PXR: ensaio Tox21 de agonismo (AID 2283988).

Reusa as funcoes do pipeline_tox21.py (leitura, rotulo, voto majoritario de replicatas) e
grava em data_pxr/, para nao regerar os 12 datasets ja curados das outras camadas.

  2283988  TOX21_PXR_LUC_Agonist  - HepG2 com luciferase sob o promotor da CYP3A4
  2283994  contra-ensaio de viabilidade do mesmo ensaio

DUAS ESCOLHAS, ambas seguindo a regra ja vigente no projeto:

1. A mascara de citotoxicidade NAO e aplicada. O projeto so a aplica aos modos
   antagonistas, porque la a morte celular imita antagonismo. Em modo agonista a perda
   de viabilidade derruba o sinal, nao o cria, entao um 'Active' agonista nao e
   explicado por citotoxicidade. As colunas citotoxico e citotox_confunde sao calculadas
   e gravadas assim mesmo, para quem quiser filtrar depois - elas so nao mexem no rotulo.
2. O contra-ensaio de autofluorescencia nao se aplica: ele e dos ensaios BLA, e este e LUC.

Uso: python pipeline_pxr_tox21.py [diretorio/com/os/csv.gz]
"""
import os
import sys

sys.path.insert(0, os.getcwd())
import numpy as np
import pandas as pd

import curation
import pipeline_tox21 as P21

AID, AID_VIAB = 2283988, 2283994
RECEPTOR, MODO, LINHAGEM = "PXR", "agonista", "HepG2-LUC"
SAIDA = "data_pxr"


# Read one PubChem assay CSV from an arbitrary directory, reusing pipeline_tox21's parser.
def ler_de(diretorio, aid):
    original = P21.DIR
    try:
        P21.DIR = diretorio
        return P21.ler(aid)
    finally:
        P21.DIR = original


if __name__ == "__main__":
    cand = [sys.argv[1]] if len(sys.argv) > 1 else []
    cand += [P21.DIR, os.path.expanduser("~/Downloads/dados/tox21")]
    diretorio = next((c for c in cand if c and os.path.exists(f"{c}/{AID}.csv.gz")), None)
    if diretorio is None:
        raise SystemExit(f"{AID}.csv.gz nao encontrado. Procurei em: {', '.join(map(str, cand))}")
    print(f"lendo de: {diretorio}", flush=True)
    os.makedirs(SAIDA, exist_ok=True)

    d = ler_de(diretorio, AID)
    d["aid"] = AID
    d["ensaio"] = "TOX21_PXR_LUC_Agonist"
    d["receptor"] = RECEPTOR
    d["modo"] = MODO
    d["linhagem"] = LINHAGEM
    d["pAC50"] = np.where(d.outcome == "Active", 6.0 - np.log10(d.AC50.where(d.AC50 > 0)), np.nan)
    d["ativo"] = (d.outcome == "Active").astype(float)
    d.loc[d.outcome == "Inconclusive", "ativo"] = np.nan
    print(f"testados {len(d)} | ativos {int((d.outcome == 'Active').sum())} "
          f"| inconclusivos {int((d.outcome == 'Inconclusive').sum())}", flush=True)

    # Viabilidade: registrada, nunca aplicada ao rotulo (ver docstring).
    if os.path.exists(f"{diretorio}/{AID_VIAB}.csv.gz"):
        v = ler_de(diretorio, AID_VIAB).rename(columns={"AC50": "AC50_viab",
                                                        "outcome": "outcome_viab"})
        d = d.merge(v[["cid", "AC50_viab", "outcome_viab"]].dropna(subset=["cid"])
                    .drop_duplicates("cid"), on="cid", how="left")
        d["citotoxico"] = d.outcome_viab == "Active"
        d["citotox_confunde"] = (d.citotoxico & d.AC50_viab.notna() & d.AC50.notna()
                                 & (d.AC50_viab <= 3.0 * d.AC50) & (d.outcome == "Active"))
        d["aid_viabilidade"] = AID_VIAB
        print(f"citotoxicos {int(d.citotoxico.sum())} | dos quais coincidem com o ativo "
              f"{int(d.citotox_confunde.sum())} (registrado, nao mascarado)", flush=True)
    else:
        d["citotoxico"] = np.nan
        d["citotox_confunde"] = np.nan
        d["aid_viabilidade"] = np.nan
    d["autofluor"] = np.nan          # contra-ensaio de autofluorescencia e dos ensaios BLA

    sp = "cache/pxr_tox21_std_lookup.pkl"
    if os.path.exists(sp):
        std = pd.read_pickle(sp)
    else:
        std = curation.standardize_frame(d.smiles, n_jobs=10)
        std.to_pickle(sp)
    std = std.rename(columns={"canonical_smiles": "smiles"})
    d = d.merge(std, on="smiles", how="left")
    print(f"SMILES unicos {std.smiles.nunique()} | rejeitados "
          f"{int(std.reject_reason.notna().sum())}", flush=True)

    d["ativo_especifico"] = d.ativo
    g0 = d[d.std_smiles.notna() & d.ativo_especifico.notna()].copy()

    # Replicatas (varios SIDs do mesmo composto) por voto majoritario; empate descartado.
    vot = g0.groupby("inchikey")["ativo_especifico"].agg(["size", "mean"])
    vot.columns = ["n_replicatas", "frac_replicatas_ativas"]
    vot["ativo_final"] = np.where(vot.frac_replicatas_ativas > 0.5, 1.0,
                         np.where(vot.frac_replicatas_ativas < 0.5, 0.0, np.nan))
    pac = g0[g0.ativo_especifico == 1].groupby("inchikey")["pAC50"].median()
    g = (g0.sort_values("ativo_especifico", ascending=False).drop_duplicates("inchikey")
           .set_index("inchikey").join(vot).join(pac.rename("pAC50_mediana")).reset_index())
    g = g[g.ativo_final.notna()].copy()
    g["ativo_especifico"] = g.ativo_final
    g["pAC50"] = np.where(g.ativo_final == 1, g.pAC50_mediana, np.nan)
    ep = f"TOX21_{RECEPTOR}_{MODO}_{LINHAGEM.split('-')[0]}"
    g["endpoint"] = ep
    g.to_csv(f"{SAIDA}/dataset_{ep}.csv", index=False)

    resumo = pd.DataFrame([{
        "aid": AID, "ensaio": "TOX21_PXR_LUC_Agonist", "receptor": RECEPTOR, "modo": MODO,
        "linhagem": LINHAGEM, "n_testados": len(d),
        "n_ativos": int((d.outcome == "Active").sum()),
        "n_inconclusivos": int((d.outcome == "Inconclusive").sum()),
        "n_citotox_confunde_registrado": int(d.citotox_confunde.fillna(False).sum()),
        "mascara_citotox_aplicada": False,
        "contra_ensaio_viabilidade": AID_VIAB,
        "compostos_validos": len(g),
        "frac_ativos": round(float(g.ativo_especifico.mean()), 4),
        "pAC50_mediana_dos_ativos": round(float(g.loc[g.ativo_especifico == 1,
                                                      "pAC50"].median()), 2)}])
    resumo.to_csv(f"{SAIDA}/resumo_pxr_tox21.csv", index=False)
    print(f"\ndataset_{ep}.csv gravado em {SAIDA}/")
    print(resumo.T.to_string(header=False))
