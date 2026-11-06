"""Tox21 layer: nuclear receptor qHTS, with viability counter-assays and
autofluorescence controls.

AC50 is in µM; pAC50 = 6 - log10(AC50_µM). AC50 is reported for all
substances, including inactive ones — it is only used when the outcome is 'Active'.
"""
import os, sys, glob
sys.path.append(os.getcwd())
import numpy as np
import pandas as pd
import curation

DIR = "tox21/2283001_2284000"

# target assay -> (receptor, mode, cell line, viability counter-assay)
ASSAYS = {
    2283418: ("ESR1", "agonist",    "HEK293T-BLA", None),
    2283419: ("ESR1", "antagonist", "HEK293T-BLA", 2283423),
    2283424: ("ESR1", "agonist",    "VM7-LUC",     None),
    2283429: ("ESR1", "antagonist", "VM7-LUC",     2283425),
    2283390: ("AR",   "agonist",    "HEK293T-BLA", None),
    2283395: ("AR",   "antagonist", "HEK293T-BLA", 2283392),
    2283393: ("AR",   "agonist",    "MDA-kb2-LUC", None),
    2283391: ("AR",   "antagonist", "MDA-kb2-LUC", 2283396),
    2283400: ("PPARG", "agonist",   "HEK293-BLA",  None),
    2283559: ("PPARG", "antagonist", "HEK293-BLA", 2283561),
    2283404: ("TR",   "agonist",    "GH3-LUC",     None),
    2283403: ("TR",   "antagonist", "GH3-LUC",     2283444),
}
AUTOFLUOR_HEK = [2283408, 2283409, 2283410]   # HEK293 blue/green/red
NAMES = {2283418: "TOX21_ERa_BLA_Agonist_ratio", 2283419: "TOX21_ERa_BLA_Antagonist_ratio",
         2283424: "TOX21_ERa_LUC_VM7_Agonist", 2283429: "TOX21_ERa_LUC_VM7_Antagonist_0.5nM_E2",
         2283390: "TOX21_AR_BLA_Agonist_ratio", 2283395: "TOX21_AR_BLA_Antagonist_ratio",
         2283393: "TOX21_AR_LUC_MDAKB2_Agonist", 2283391: "TOX21_AR_LUC_MDAKB2_Antagonist_10nM_R1881",
         2283400: "TOX21_PPARg_BLA_Agonist_ratio", 2283559: "TOX21_PPARg_BLA_Antagonist_ratio",
         2283404: "TOX21_TR_LUC_GH3_Agonist", 2283403: "TOX21_TR_LUC_GH3_Antagonist"}


# Read one PubChem assay CSV, dropping the two metadata header rows and coercing the numeric columns.
def read(aid):
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
            d = read(a)
            af |= set(d.loc[d.outcome == "Active", "cid"].dropna().astype(int))
    print("substances with autofluorescence in HEK293:", len(af))

    summary, all_items = [], []
    for aid, (rec, mode, cel, viab) in ASSAYS.items():
        d = read(aid)
        d["aid"] = aid; d["assay"] = NAMES[aid]; d["receptor"] = rec
        d["mode"] = mode; d["cell_line"] = cel
        d["pAC50"] = np.where(d.outcome == "Active", 6.0 - np.log10(d.AC50.where(d.AC50 > 0)), np.nan)
        d["active"] = (d.outcome == "Active").astype(float)
        d.loc[d.outcome == "Inconclusive", "active"] = np.nan
        if viab is not None and os.path.exists(f"{DIR}/{viab}.csv.gz"):
            v = read(viab).rename(columns={"AC50": "AC50_viability", "outcome": "outcome_viability"})
            d = d.merge(v[["cid", "AC50_viability", "outcome_viability"]].dropna(subset=["cid"]).drop_duplicates("cid"),
                        on="cid", how="left")
            d["cytotoxic"] = (d.outcome_viability == "Active")
            d["cytotox_confounds"] = (d.cytotoxic & d.AC50_viability.notna() & d.AC50.notna()
                                     & (d.AC50_viability <= 3.0 * d.AC50) & (d.outcome == "Active"))
            d["aid_viability"] = viab
        else:
            d["cytotoxic"] = np.nan; d["cytotox_confounds"] = np.nan; d["aid_viability"] = np.nan
        d["autofluorescence"] = d.cid.isin(af) if cel.endswith("BLA") else np.nan
        all_items.append(d)
        summary.append({"aid": aid, "assay": NAMES[aid], "receptor": rec, "mode": mode,
                       "cell_line": cel, "n_tested": len(d),
                       "n_active": int((d.outcome == "Active").sum()),
                       "n_inconclusive": int((d.outcome == "Inconclusive").sum()),
                       "n_cytotox_confounds": int(d.cytotox_confounds.fillna(False).sum()),
                       "n_autofluorescence_active": int((d.autofluorescence.fillna(False) & (d.outcome == "Active")).sum()),
                       "viability_counter_assay": viab})
    tox = pd.concat(all_items, ignore_index=True)

    sp = "cache/tox21_std_lookup.pkl"
    if os.path.exists(sp):
        std = pd.read_pickle(sp)
    else:
        std = curation.standardize_frame(tox.smiles, n_jobs=10)
        std.to_pickle(sp)
    std = std.rename(columns={"canonical_smiles": "smiles"})
    tox = tox.merge(std, on="smiles", how="left")
    print("Unique SMILES:", std.smiles.nunique(), "| rejected:", int(std.reject_reason.notna().sum()))
    print(std.reject_reason.str.split(":").str[0].value_counts().head(8).to_dict())

    # Training label: an "active" whose activity coincides with loss of viability
    # is not distinguishable from cell death — becomes NaN (stays in the CSV with the flag), not inactive.
    conf = (tox.cytotox_confounds.fillna(False).astype(bool)
            | (tox.autofluorescence.fillna(False).astype(bool) & (tox.active == 1)))
    tox["specific_active"] = tox.active.where(~((tox.active == 1) & conf))
    tox.to_pickle("cache/tox21_curated.pkl")

    # Tox21 tests several samples (SIDs) of the same compound. Replicates are
    # resolved by majority vote among the conclusive ones; ties are discarded.
    for aid, (rec, mode, cel, viab) in ASSAYS.items():
        g0 = tox[(tox.aid == aid) & tox.std_smiles.notna() & tox.specific_active.notna()].copy()
        n_voting_models = g0.groupby("inchikey")["specific_active"].agg(["size", "mean"])
        n_voting_models.columns = ["n_replicates", "active_replicate_fraction"]
        n_voting_models["final_active"] = np.where(n_voting_models.active_replicate_fraction > 0.5, 1.0,
                             np.where(n_voting_models.active_replicate_fraction < 0.5, 0.0, np.nan))
        pac = g0[g0.specific_active == 1].groupby("inchikey")["pAC50"].median()
        g = (g0.sort_values("specific_active", ascending=False).drop_duplicates("inchikey")
               .set_index("inchikey").join(n_voting_models).join(pac.rename("pAC50_median")).reset_index())
        g = g[g.final_active.notna()].copy()
        g["specific_active"] = g["final_active"]
        g["pAC50"] = np.where(g.final_active == 1, g.pAC50_median, np.nan)
        ep = f"TOX21_{rec}_{mode}_{cel.split('-')[0]}"
        g["endpoint"] = ep
        g.to_csv(f"data_tox21/dataset_{ep}.csv", index=False)
        for r in summary:
            if r["aid"] == aid:
                r["compounds_valid"] = len(g)
                r["frac_actives"] = round(float(g.specific_active.mean()), 4)
    rs = pd.DataFrame(summary)
    rs.to_csv("data_tox21/tox21_summary.csv", index=False)
    print("\n", rs.to_string(index=False))


if __name__ == "__main__":
    main()
