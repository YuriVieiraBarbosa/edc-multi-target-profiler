"""Camada Tox21: qHTS de receptores nucleares, com contra-ensaios de
viabilidade e controles de autofluorescência.

AC50 vem em µM; pAC50 = 6 - log10(AC50_µM). A AC50 é reportada para todas as
substâncias, inclusive inativas — só é usada quando o outcome é 'Active'.
"""
import os, sys, glob
sys.path.append(os.getcwd())
import numpy as np
import pandas as pd
import curation

DIR = "tox21/2283001_2284000"

# ensaio-alvo -> (receptor, modo, linhagem, contra-ensaio de viabilidade)
ASSAYS = {
    2283418: ("ESR1", "agonista",    "HEK293T-BLA", None),
    2283419: ("ESR1", "antagonista", "HEK293T-BLA", 2283423),
    2283424: ("ESR1", "agonista",    "VM7-LUC",     None),
    2283429: ("ESR1", "antagonista", "VM7-LUC",     2283425),
    2283390: ("AR",   "agonista",    "HEK293T-BLA", None),
    2283395: ("AR",   "antagonista", "HEK293T-BLA", 2283392),
    2283393: ("AR",   "agonista",    "MDA-kb2-LUC", None),
    2283391: ("AR",   "antagonista", "MDA-kb2-LUC", 2283396),
    2283400: ("PPARG", "agonista",   "HEK293-BLA",  None),
    2283559: ("PPARG", "antagonista", "HEK293-BLA", 2283561),
    2283404: ("TR",   "agonista",    "GH3-LUC",     None),
    2283403: ("TR",   "antagonista", "GH3-LUC",     2283444),
}
AUTOFLUOR_HEK = [2283408, 2283409, 2283410]   # HEK293 azul/verde/vermelho
NOMES = {2283418: "TOX21_ERa_BLA_Agonist_ratio", 2283419: "TOX21_ERa_BLA_Antagonist_ratio",
         2283424: "TOX21_ERa_LUC_VM7_Agonist", 2283429: "TOX21_ERa_LUC_VM7_Antagonist_0.5nM_E2",
         2283390: "TOX21_AR_BLA_Agonist_ratio", 2283395: "TOX21_AR_BLA_Antagonist_ratio",
         2283393: "TOX21_AR_LUC_MDAKB2_Agonist", 2283391: "TOX21_AR_LUC_MDAKB2_Antagonist_10nM_R1881",
         2283400: "TOX21_PPARg_BLA_Agonist_ratio", 2283559: "TOX21_PPARg_BLA_Antagonist_ratio",
         2283404: "TOX21_TR_LUC_GH3_Agonist", 2283403: "TOX21_TR_LUC_GH3_Antagonist"}


# Read one PubChem assay CSV, dropping the two metadata header rows and coercing the numeric columns.
def ler(aid):
    d = pd.read_csv(f"{DIR}/{aid}.csv.gz", low_memory=False)
    d = d[pd.to_numeric(d.PUBCHEM_SID, errors="coerce").notna()].copy()
    for c in ("AC50", "HITC", "PUBCHEM_CID"):
        if c in d:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.rename(columns={"PUBCHEM_CID": "cid", "PUBCHEM_SID": "sid",
                          "PUBCHEM_EXT_DATASOURCE_SMILES": "smiles",
                          "PUBCHEM_ACTIVITY_OUTCOME": "outcome"})
    return d[["sid", "cid", "smiles", "outcome", "AC50", "HITC"]]


# Build one curated dataset per Tox21 assay, applying the cytotoxicity and autofluorescence counter-screens.
def main():
    os.makedirs("data_tox21", exist_ok=True)
    af = set()
    for a in AUTOFLUOR_HEK:
        if os.path.exists(f"{DIR}/{a}.csv.gz"):
            d = ler(a)
            af |= set(d.loc[d.outcome == "Active", "cid"].dropna().astype(int))
    print("substâncias com autofluorescência em HEK293:", len(af))

    resumo, todos = [], []
    for aid, (rec, modo, cel, viab) in ASSAYS.items():
        d = ler(aid)
        d["aid"] = aid; d["ensaio"] = NOMES[aid]; d["receptor"] = rec
        d["modo"] = modo; d["linhagem"] = cel
        d["pAC50"] = np.where(d.outcome == "Active", 6.0 - np.log10(d.AC50.where(d.AC50 > 0)), np.nan)
        d["ativo"] = (d.outcome == "Active").astype(float)
        d.loc[d.outcome == "Inconclusive", "ativo"] = np.nan
        if viab is not None and os.path.exists(f"{DIR}/{viab}.csv.gz"):
            v = ler(viab).rename(columns={"AC50": "AC50_viab", "outcome": "outcome_viab"})
            d = d.merge(v[["cid", "AC50_viab", "outcome_viab"]].dropna(subset=["cid"]).drop_duplicates("cid"),
                        on="cid", how="left")
            d["citotoxico"] = (d.outcome_viab == "Active")
            d["citotox_confunde"] = (d.citotoxico & d.AC50_viab.notna() & d.AC50.notna()
                                     & (d.AC50_viab <= 3.0 * d.AC50) & (d.outcome == "Active"))
            d["aid_viabilidade"] = viab
        else:
            d["citotoxico"] = np.nan; d["citotox_confunde"] = np.nan; d["aid_viabilidade"] = np.nan
        d["autofluor"] = d.cid.isin(af) if cel.endswith("BLA") else np.nan
        todos.append(d)
        resumo.append({"aid": aid, "ensaio": NOMES[aid], "receptor": rec, "modo": modo,
                       "linhagem": cel, "n_testados": len(d),
                       "n_ativos": int((d.outcome == "Active").sum()),
                       "n_inconclusivos": int((d.outcome == "Inconclusive").sum()),
                       "n_citotox_confunde": int(d.citotox_confunde.fillna(False).sum()),
                       "n_autofluor_ativo": int((d.autofluor.fillna(False) & (d.outcome == "Active")).sum()),
                       "contra_ensaio_viabilidade": viab})
    tox = pd.concat(todos, ignore_index=True)

    sp = "cache/tox21_std_lookup.pkl"
    if os.path.exists(sp):
        std = pd.read_pickle(sp)
    else:
        std = curation.standardize_frame(tox.smiles, n_jobs=10)
        std.to_pickle(sp)
    std = std.rename(columns={"canonical_smiles": "smiles"})
    tox = tox.merge(std, on="smiles", how="left")
    print("SMILES únicos:", std.smiles.nunique(), "| rejeitados:", int(std.reject_reason.notna().sum()))
    print(std.reject_reason.str.split(":").str[0].value_counts().head(8).to_dict())

    # Rótulo de treino: um "ativo" cuja atividade coincide com perda de viabilidade
    # não é distinguível de morte celular — vira NaN (fica no CSV com a flag), não inativo.
    conf = (tox.citotox_confunde.fillna(False).astype(bool)
            | (tox.autofluor.fillna(False).astype(bool) & (tox.ativo == 1)))
    tox["ativo_especifico"] = tox.ativo.where(~((tox.ativo == 1) & conf))
    tox.to_pickle("cache/tox21_curado.pkl")

    # O Tox21 testa várias amostras (SIDs) do mesmo composto. Réplicas são
    # resolvidas por voto majoritário entre as conclusivas; empate é descartado.
    for aid, (rec, modo, cel, viab) in ASSAYS.items():
        g0 = tox[(tox.aid == aid) & tox.std_smiles.notna() & tox.ativo_especifico.notna()].copy()
        vot = g0.groupby("inchikey")["ativo_especifico"].agg(["size", "mean"])
        vot.columns = ["n_replicatas", "frac_replicatas_ativas"]
        vot["ativo_final"] = np.where(vot.frac_replicatas_ativas > 0.5, 1.0,
                             np.where(vot.frac_replicatas_ativas < 0.5, 0.0, np.nan))
        pac = g0[g0.ativo_especifico == 1].groupby("inchikey")["pAC50"].median()
        g = (g0.sort_values("ativo_especifico", ascending=False).drop_duplicates("inchikey")
               .set_index("inchikey").join(vot).join(pac.rename("pAC50_mediana")).reset_index())
        g = g[g.ativo_final.notna()].copy()
        g["ativo_especifico"] = g["ativo_final"]
        g["pAC50"] = np.where(g.ativo_final == 1, g.pAC50_mediana, np.nan)
        ep = f"TOX21_{rec}_{modo}_{cel.split('-')[0]}"
        g["endpoint"] = ep
        g.to_csv(f"data_tox21/dataset_{ep}.csv", index=False)
        for r in resumo:
            if r["aid"] == aid:
                r["compostos_validos"] = len(g)
                r["frac_ativos"] = round(float(g.ativo_especifico.mean()), 4)
    rs = pd.DataFrame(resumo)
    rs.to_csv("data_tox21/resumo_tox21.csv", index=False)
    print("\n", rs.to_string(index=False))


if __name__ == "__main__":
    main()
