"""PXR screening axis: Tox21 agonism assay (AID 2283988).

Reuses the functions from pipeline_tox21.py (reading, label, majority vote across replicates) and
writes to data_pxr/, so as not to regenerate the 12 datasets already curated for the other layers.

  2283988  TOX21_PXR_LUC_Agonist  - HepG2 with luciferase under the CYP3A4 promoter
  2283994  viability counter-assay for the same assay

TWO CHOICES, both following the rule already in force in the project:

1. The cytotoxicity mask is NOT applied. The project only applies it to antagonist
   modes, because there cell death mimics antagonism. In agonist mode, loss of
   viability suppresses the signal rather than creating it, so an agonist 'Active' is not
   explained by cytotoxicity. The cytotoxic and cytotox_confounds columns are calculated
   and saved anyway, for anyone who wants to filter afterward - they just don't affect the label.
2. The autofluorescence counter-assay does not apply: it belongs to the BLA assays, and this is LUC.

Usage: python pipeline_pxr_tox21.py [directory/with/the/csv.gz]
"""
import os
import sys

sys.path.insert(0, os.getcwd())
import numpy as np
import pandas as pd

import curation
import pipeline_tox21 as P21

AID, AID_VIABILITY = 2283988, 2283994
RECEPTOR, MODE, CELL_LINE = "PXR", "agonist", "HepG2-LUC"
OUTPUT_DIR = "data_pxr"


# Read one PubChem assay CSV from an arbitrary directory, reusing pipeline_tox21's parser.
def read_from(directory, aid):
    original = P21.DIR
    try:
        P21.DIR = directory
        return P21.read(aid)
    finally:
        P21.DIR = original


if __name__ == "__main__":
    cand = [sys.argv[1]] if len(sys.argv) > 1 else []
    cand += [P21.DIR, os.path.expanduser("~/Downloads/data/tox21")]
    directory = next((c for c in cand if c and os.path.exists(f"{c}/{AID}.csv.gz")), None)
    if directory is None:
        raise SystemExit(f"{AID}.csv.gz not found. Looked in: {', '.join(map(str, cand))}")
    print(f"reading from: {directory}", flush=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    d = read_from(directory, AID)
    d["aid"] = AID
    d["assay"] = "TOX21_PXR_LUC_Agonist"
    d["receptor"] = RECEPTOR
    d["mode"] = MODE
    d["cell_line"] = CELL_LINE
    d["pAC50"] = np.where(d.outcome == "Active", 6.0 - np.log10(d.AC50.where(d.AC50 > 0)), np.nan)
    d["active"] = (d.outcome == "Active").astype(float)
    d.loc[d.outcome == "Inconclusive", "active"] = np.nan
    print(f"tested {len(d)} | actives {int((d.outcome == 'Active').sum())} "
          f"| inconclusive {int((d.outcome == 'Inconclusive').sum())}", flush=True)

    # Viability: recorded, never applied to the label (see docstring).
    if os.path.exists(f"{directory}/{AID_VIABILITY}.csv.gz"):
        v = read_from(directory, AID_VIABILITY).rename(columns={"AC50": "AC50_viability",
                                                        "outcome": "outcome_viability"})
        d = d.merge(v[["cid", "AC50_viability", "outcome_viability"]].dropna(subset=["cid"])
                    .drop_duplicates("cid"), on="cid", how="left")
        d["cytotoxic"] = d.outcome_viability == "Active"
        d["cytotox_confounds"] = (d.cytotoxic & d.AC50_viability.notna() & d.AC50.notna()
                                 & (d.AC50_viability <= 3.0 * d.AC50) & (d.outcome == "Active"))
        d["aid_viability"] = AID_VIABILITY
        print(f"cytotoxic {int(d.cytotoxic.sum())} | of which coincide with the active "
              f"{int(d.cytotox_confounds.sum())} (recorded, not masked)", flush=True)
    else:
        d["cytotoxic"] = np.nan
        d["cytotox_confounds"] = np.nan
        d["aid_viability"] = np.nan
    d["autofluorescence"] = np.nan          # autofluorescence counter-assay belongs to the BLA assays

    sp = "cache/pxr_tox21_std_lookup.pkl"
    if os.path.exists(sp):
        std = pd.read_pickle(sp)
    else:
        std = curation.standardize_frame(d.smiles, n_jobs=10)
        std.to_pickle(sp)
    std = std.rename(columns={"canonical_smiles": "smiles"})
    d = d.merge(std, on="smiles", how="left")
    print(f"Unique SMILES {std.smiles.nunique()} | rejected "
          f"{int(std.reject_reason.notna().sum())}", flush=True)

    d["specific_active"] = d.active
    g0 = d[d.std_smiles.notna() & d.specific_active.notna()].copy()

    # Replicates (multiple SIDs of the same compound) resolved by majority vote; ties discarded.
    n_voting_models = g0.groupby("inchikey")["specific_active"].agg(["size", "mean"])
    n_voting_models.columns = ["n_replicates", "active_replicate_fraction"]
    n_voting_models["final_active"] = np.where(n_voting_models.active_replicate_fraction > 0.5, 1.0,
                         np.where(n_voting_models.active_replicate_fraction < 0.5, 0.0, np.nan))
    pac = g0[g0.specific_active == 1].groupby("inchikey")["pAC50"].median()
    g = (g0.sort_values("specific_active", ascending=False).drop_duplicates("inchikey")
           .set_index("inchikey").join(n_voting_models).join(pac.rename("pAC50_median")).reset_index())
    g = g[g.final_active.notna()].copy()
    g["specific_active"] = g.final_active
    g["pAC50"] = np.where(g.final_active == 1, g.pAC50_median, np.nan)
    ep = f"TOX21_{RECEPTOR}_{MODE}_{CELL_LINE.split('-')[0]}"
    g["endpoint"] = ep
    g.to_csv(f"{OUTPUT_DIR}/dataset_{ep}.csv", index=False)

    summary = pd.DataFrame([{
        "aid": AID, "assay": "TOX21_PXR_LUC_Agonist", "receptor": RECEPTOR, "mode": MODE,
        "cell_line": CELL_LINE, "n_tested": len(d),
        "n_active": int((d.outcome == "Active").sum()),
        "n_inconclusive": int((d.outcome == "Inconclusive").sum()),
        "n_cytotox_confounds_registered": int(d.cytotox_confounds.fillna(False).sum()),
        "applied_cytotox_mask": False,
        "viability_counter_assay": AID_VIABILITY,
        "compounds_valid": len(g),
        "frac_actives": round(float(g.specific_active.mean()), 4),
        "pAC50_median_of_actives": round(float(g.loc[g.specific_active == 1,
                                                      "pAC50"].median()), 2)}])
    summary.to_csv(f"{OUTPUT_DIR}/summary_pxr_tox21.csv", index=False)
    print(f"\ndataset_{ep}.csv saved to {OUTPUT_DIR}/")
    print(summary.T.to_string(header=False))
