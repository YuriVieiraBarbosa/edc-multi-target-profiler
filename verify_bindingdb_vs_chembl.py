"""Checks the BindingDB layer against ChEMBL on the compounds that both sources have.

A new layer is only trustworthy after it matches the already curated source on the
shared points. In this project this check already paid off: in a previous source it revealed a
group of pairs that diverged by a factor of a thousand, concentrated in a few origins, with
a unit-error signature.

What it looks for, per pair (compound, target, measure) present in both sources:
  |diff| <= 1,0    inside the replicate tolerance that the funnel itself uses
  |diff| > 1,0     real divergence, listed
  |diff| ~ 3,0     candidate for a unit error (factor 1.000, typically nM x uM)

Usage:  python verify_bindingdb_vs_chembl.py [path/to/dataset_IC50_curated.csv]
Output: data_bindingdb/check_vs_chembl.csv and the summary on screen.

The ChEMBL IC50 layer is not versioned in the folder; pass its path as
an argument to include it, otherwise the check covers only Ki, Kd and EC50.
"""
import glob
import os
import sys

sys.path.insert(0, os.getcwd())
import numpy as np
import pandas as pd
from scipy import stats

OUTPUT_DIR = "data_bindingdb"


# Load the curated ChEMBL side: the Ki/Kd/EC50 layers plus, if given, the IC50 dataset.
def side_chembl(ic50_path=None):
    parts = []
    folders = ["datasets_layers_by_endpoint", "layers_curated_by_endpoint", "data_layers"]
    for p in folders:
        for a in glob.glob(os.path.join(p, "dataset_*.csv")):
            d = pd.read_csv(a, usecols=["inchikey", "endpoint", "pAct_final"])
            d["target"] = d.endpoint.str.split("_").str[0]
            d["measure"] = d.endpoint.str.split("_").str[1]
            parts.append(d[["inchikey", "target", "measure", "pAct_final"]])
        if parts:
            print(f"ChEMBL layers read from: {p}")
            break
    if ic50_path and os.path.exists(ic50_path):
        d = pd.read_csv(ic50_path, low_memory=False, usecols=["inchikey", "endpoint", "pIC50_final"])
        d["target"] = d.endpoint.str.split("_").str[0]
        d["measure"] = "IC50"
        parts.append(d.rename(columns={"pIC50_final": "pAct_final"})[
            ["inchikey", "target", "measure", "pAct_final"]])
        print(f"IC50 layer read from: {ic50_path}")
    else:
        print("WARNING: without the ChEMBL IC50 dataset; the check covers only Ki, Kd and EC50")
    if not parts:
        raise SystemExit("no curated ChEMBL dataset found")
    return pd.concat(parts, ignore_index=True)


# Load the curated BindingDB side written by pipeline_bindingdb.py.
def side_bindingdb():
    fnames = sorted(glob.glob(f"{OUTPUT_DIR}/dataset_*_BDB*.csv"))
    if not fnames:
        raise SystemExit(f"no dataset in {OUTPUT_DIR}/. Run pipeline_bindingdb.py first.")
    parts = []
    for a in fnames:
        d = pd.read_csv(a, low_memory=False,
                        usecols=["inchikey", "endpoint", "pAct_final", "target_label",
                                 "standard_type", "provenance", "sources"])
        parts.append(d)
    return pd.concat(parts, ignore_index=True)


if __name__ == "__main__":
    ch = side_chembl(sys.argv[1] if len(sys.argv) > 1 else None)
    bd = side_bindingdb()
    print(f"\ncurated pairs: ChEMBL {len(ch)} | BindingDB {len(bd)}")

    m = bd.merge(ch, left_on=["inchikey", "target_label", "standard_type"],
                 right_on=["inchikey", "target", "measure"],
                 suffixes=("_bdb", "_chembl"))
    m = (m.groupby(["inchikey", "target_label", "standard_type", "provenance", "sources"],
                   as_index=False)
           .agg(pAct_bdb=("pAct_final_bdb", "median"),
                pAct_chembl=("pAct_final_chembl", "median")))
    m['diff'] = m.pAct_bdb - m.pAct_chembl
    m["suspect_unit"] = m['diff'].abs().between(2.5, 3.5)
    m.to_csv(f"{OUTPUT_DIR}/check_vs_chembl.csv", index=False)

    if not len(m):
        print("\nno pair in common: nothing to check")
        sys.exit(0)

    print(f"\n=== {len(m)} pairs (compound, target, measure) in both sources ===")
    print(f"  median difference BindingDB - ChEMBL : {m['diff'].median():+.2f} log")
    print(f"  |diff| median                        : {m['diff'].abs().median():.2f} log")
    print(f"  inside the 1,0 log tolerance      : {int((m['diff'].abs() <= 1).sum())} "
          f"({(m['diff'].abs() <= 1).mean():.1%})")
    print(f"  |diff| > 1,0                          : {int((m['diff'].abs() > 1).sum())}")
    print(f"  |diff| between 2,5 and 3,5 (factor ~1.000) : {int(m.suspect_unit.sum())}")
    if len(m) > 2:
        print(f"  Spearman between the two sources         : {stats.spearmanr(m.pAct_bdb, m.pAct_chembl).statistic:.2f}")

    print("\n=== by origin ===")
    print(m.groupby("provenance").agg(
        pairs=("diff", "size"), median_diff=("diff", "median"),
        absdif_median=("diff", lambda s: s.abs().median()),
        above_1=("diff", lambda s: int((s.abs() > 1).sum())),
        suspect_unit=("suspect_unit", "sum")).round(2).to_string())

    suspects = m[m.suspect_unit]
    if len(suspects):
        print(f"\n=== {len(suspects)} pairs with unit error signature ===")
        print(suspects.nlargest(min(12, len(suspects)), "diff", keep="all")[
            ["inchikey", "target_label", "standard_type", "provenance",
             "pAct_bdb", "pAct_chembl", "diff", "sources"]].round(2).to_string(index=False))
    print(f"\nwritten: {OUTPUT_DIR}/check_vs_chembl.csv")
