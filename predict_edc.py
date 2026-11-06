#!/usr/bin/env python
"""Multi-target prediction of endocrine-disrupting potential from SMILES.

Usage:
    python predict_edc.py --smiles "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1" "OCCOC(=O)..."
    python predict_edc.py --input molecules.csv --smiles-col smiles --out predictions.csv

For each trained endpoint, returns the predicted pIC50, probability of being active
(pIC50 >= 6, i.e. IC50 <= 1 uM), the Tanimoto similarity to the nearest
training neighbor, and the applicability domain flag. Predictions outside the
applicability domain must NOT be used for prioritization.
"""
import argparse, glob, os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import joblib
import curation
import features as F

LABELS = {"AR": "Androgen receptor", "ESR1": "Estrogen receptor alpha",
           "ESR2": "Estrogen receptor beta", "PPARG": "PPAR-gamma",
           "THRA": "Thyroid hormone receptor alpha",
           "THRB": "Thyroid hormone receptor beta"}
ASSAY = {"B": "binding", "F": "functional"}


# Load the regression models of the earlier single-layer run.
def load_models(model_dir="models"):
    mods = {}
    for p in sorted(glob.glob(os.path.join(model_dir, "*_regression.joblib"))):
        ep = os.path.basename(p).replace("_regression.joblib", "")
        mods[ep] = {"reg": joblib.load(p)}
        pc = os.path.join(model_dir, f"{ep}_classification.joblib")
        if os.path.exists(pc):
            mods[ep]["clf"] = joblib.load(pc)
    return mods


# Predict pIC50 for a list of SMILES using those models.
def predict(smiles_list, model_dir="models"):
    mods = load_models(model_dir)
    if not mods:
        raise SystemExit(f"no model found in {model_dir}/")
    std = [curation.standardize(s) for s in smiles_list]
    base = pd.DataFrame({"input_smiles": smiles_list,
                         "std_smiles": [d["std_smiles"] for d in std],
                         "inchikey": [d["inchikey"] for d in std],
                         "rejection_reason": [d["reject_reason"] for d in std]})
    val = base["std_smiles"].notna().to_numpy()
    if val.sum() == 0:
        return base
    X, bits, _, ok = F.featurize(base.loc[val, "std_smiles"].tolist())
    idx = base.index[val][ok]
    for ep, m in sorted(mods.items()):
        target, kind = ep.rsplit("_", 1)
        pref = f"{target}_{kind}"
        base.loc[idx, f"{pref}_pIC50_prev"] = m["reg"]["model"].predict(X).round(2)
        sim = F.tanimoto_max_to_train(bits, m["reg"]["train_bits"])
        base.loc[idx, f"{pref}_tanimoto_nn"] = sim.round(3)
        base.loc[idx, f"{pref}_in_domain"] = sim >= m["reg"]["ad_threshold"]
        if "clf" in m:
            base.loc[idx, f"{pref}_prob_active"] = m["clf"]["model"].predict_proba(X)[:, 1].round(3)
    probs = [c for c in base.columns if c.endswith("_prob_active")]
    domain_cols = [c.replace("_prob_active", "_in_domain") for c in probs]
    active_in_domain = np.zeros(len(base))
    for pc, dc in zip(probs, domain_cols):
        active_in_domain += ((base[pc] >= 0.5) & (base[dc] == True)).astype(int)
    base["n_endpoints_active_in_domain"] = active_in_domain.astype(int)
    base["multitarget_profile"] = np.where(base["std_smiles"].isna(), "rejected structure",
                                        np.where(active_in_domain >= 3, "multi-target (>=3 endpoints)",
                                        np.where(active_in_domain >= 1, "selective (1-2 endpoints)",
                                                 "no alert")))
    return base


# Command-line entry point of the first predictor, kept for reference only.
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smiles", nargs="+", help="one or more SMILES")
    ap.add_argument("--input", help="Input CSV")
    ap.add_argument("--smiles-col", default="smiles")
    ap.add_argument("--name-col", default=None)
    ap.add_argument("--model-dir", default="models")
    ap.add_argument("--out", default="predictions_edc.csv")
    a = ap.parse_args()
    if a.input:
        df = pd.read_csv(a.input)
        smi = df[a.smiles_col].astype(str).tolist()
        names = df[a.name_col].tolist() if a.name_col else None
    elif a.smiles:
        smi, names = a.smiles, None
    else:
        ap.error("provide --smiles or --input")
    out = predict(smi, a.model_dir)
    if names is not None:
        out.insert(0, "name", names)
    out.to_csv(a.out, index=False)
    print(f"{len(out)} molecules -> {a.out}")
    cols = ["name"] if names is not None else ["input_smiles"]
    cols += ["n_endpoints_active_in_domain", "multitarget_profile"]
    print(out[[c for c in cols if c in out]].to_string(index=False))


if __name__ == "__main__":
    main()
