"""ChEMBL layer for the pregnane X receptor (PXR / NR1I2, UniProt O75469).

PXR enters the project as a xenosensor: it regulates CYP3A4 and mediates much of
the effect of endocrine disruptors through enzymatic induction, besides sharing ligands with
ER, AR, and PPARgamma. Same funnel, same standardization, and same endpoint grammar as the
other ChEMBL layers.

ONE DESIGN DIFFERENCE, forced by what the source provides. In PXR, ChEMBL's `assay_type` field
poorly separates binding from functional - a defect already noted in the project for
other targets, but here it is extreme:

  kind A (ADMET)  421 records, the LARGEST group: per the description these are CYP3A4 reporter assays in
                  HepG2/DPX-2, i.e. FUNCTIONAL, plus one binding TR-FRET
  kind B          646 records: mixes BODIPY FL vindoline displacement (real
                  binding) with a luciferase reporter in HepG2 (functional)
  kind F           66 records, reporter
  kind T           32 records, transactivation

Discarding kind A, as the funnel does for the other layers, would throw away 36% of the records
of the target. Merging it into B would violate the project rule of never mixing
incompatible assays. The adopted output is the third option: each assay_type becomes its OWN endpoint
(PXR_EC50_A, PXR_IC50_B, ...), which is the same rule already applied to B and F, just extended.
The consequence must travel with any reading: no PXR endpoint can be
called 'binding' or 'functional' by its suffix.

Usage:  python pipeline_pxr.py [path/to/chembl_37.db]
Output: data_pxr/dataset_PXR_<MEASURE>_<KIND>.csv, funnel log, and summary per endpoint.
"""
import os
import re
import sqlite3
import sys

sys.path.insert(0, os.getcwd())
import numpy as np
import pandas as pd

import curation
from extract_chembl import ACT_SQL, TARGET_SQL

ACC_PXR = "O75469"
LABEL = "PXR"
TYPES = ["IC50", "Ki", "Kd", "EC50"]
MIN_REG = 50
OUTPUT_DIR = "data_pxr"

BAD_VALIDITY = {"Outside typical range", "Non standard unit for type",
                "Potential missing data", "Potential transcription error",
                "Potential author error", "Author confirmed error", "Potential duplicate"}
UNIT_TO_NM = {"nM": 1.0, "uM": 1e3, "µM": 1e3, "mM": 1e6, "M": 1e9, "pM": 1e-3}
AGONIST_RE = re.compile(r"\bagonis", re.I)
ANTAG_RE = re.compile(r"\bantagonis|\binhibit|\bblock", re.I)
PROV = ["molecule_chembl_id", "molecule_name", "canonical_smiles", "std_smiles",
        "inchikey", "standard_type", "standard_relation", "standard_value",
        "standard_units", "value_nM", "pAct", "pchembl_value", "assay_chembl_id",
        "assay_description", "assay_type", "assay_mode", "assay_organism",
        "confidence_score", "bao_format", "target_chembl_id", "target_name",
        "accession", "target_label", "target_organism", "document_chembl_id",
        "doi", "pubmed_id", "year", "journal", "activity_source", "undefined_stereo", "mw"]


# Locate the ChEMBL SQLite dump: explicit argument, then the usual places.
def find_db(argument=None):
    cand = [argument] if argument else []
    cand += ["chembl/chembl_37/chembl_37_sqlite/chembl_37.db",
             os.path.expanduser("~/Downloads/data/chembl_37/chembl_37/"
                                "chembl_37_sqlite/chembl_37.db")]
    for c in cand:
        if c and os.path.exists(c):
            return c
    raise SystemExit("chembl_37.db not found. Pass the path as an argument.\n"
                     "Looked in: " + ", ".join(str(c) for c in cand if c))


if __name__ == "__main__":
    db = find_db(sys.argv[1] if len(sys.argv) > 1 else None)
    print(f"reading: {db}", flush=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs("cache", exist_ok=True)

    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    tg = pd.read_sql(TARGET_SQL.format(placeholders="?"), con, params=[ACC_PXR])
    single = tg[(tg.target_type == "SINGLE PROTEIN") & (tg.organism == "Homo sapiens")]
    if not len(single):
        raise SystemExit(f"no human SINGLE PROTEIN target for {ACC_PXR}")
    tids = sorted(set(single.tid))
    print(f"target: {single.target_chembl_id.iloc[0]} "
          f"({single.target_name.iloc[0]}) | tid {tids}", flush=True)

    sql = ACT_SQL.format(tids=",".join("?" * len(tids))).replace(
        "AND act.standard_type = 'IC50'",
        "AND act.standard_type IN ({})".format(",".join("?" * len(TYPES))))
    raw = pd.read_sql(sql, con, params=tids + TYPES)
    con.close()
    raw = raw.merge(single[["target_chembl_id", "accession"]].drop_duplicates(),
                    on="target_chembl_id", how="left")
    raw["target_label"] = LABEL

    log = []

    # Record the record count after one filtering step, for the curation funnel log.
    def note(step, df):
        log.append({"step": step, "records": len(df)})
        print(f"{step:50s} n={len(df):7d}", flush=True)

    note("0. raw (IC50/Ki/Kd/EC50, human PXR)", raw)
    d = raw[raw.target_organism == "Homo sapiens"]
    d = d[d.confidence_score >= 8]
    note("1. confidence_score >= 8", d)
    d = d[d.standard_relation == "="]
    note("2. relation '='", d)
    d = d[d.standard_units.isin(UNIT_TO_NM)].copy()
    d["value_nM"] = pd.to_numeric(d.standard_value, errors="coerce") * d.standard_units.map(UNIT_TO_NM)
    note("3. unit convertible to nM", d)
    d = d[d.value_nM.between(1e-3, 1e7)]
    note("4. plausible band", d)
    d = d[~d.data_validity_comment.isin(BAD_VALIDITY)]
    note("5. validity ok", d)
    d = d[d.canonical_smiles.notna() & (d.structure_type == "MOL")]
    note("6. with structure", d)

    std_path = "cache/pxr_std_lookup.pkl"
    if os.path.exists(std_path):
        std = pd.read_pickle(std_path)
    else:
        std = curation.standardize_frame(d.canonical_smiles, n_jobs=10)
        std.to_pickle(std_path)
    d = d.merge(std, on="canonical_smiles", how="left")
    d = d[d.reject_reason.isna()]
    note("7. standardized structure valid", d)

    d["pAct"] = curation.to_pic50(d.value_nM)
    # The mode comes from the description, as in the other layers. The regex order tests antagonist
    # first, then 'inhibit'/'block' win over 'agonis' - a known limitation.
    desc = d.assay_description.fillna("")
    d["assay_mode"] = np.where(desc.str.contains(ANTAG_RE), "antagonist",
                      np.where(desc.str.contains(AGONIST_RE), "agonist", "undefined"))
    d["source"] = "ChEMBL37"
    d.to_pickle("cache/pxr_filtered.pkl")

    summary = []
    for (kind, atype), grp in d.groupby(["standard_type", "assay_type"]):
        if len(grp) < MIN_REG:
            continue
        agg = curation.aggregate_replicates(grp, ["inchikey"], "pAct", max_spread=1.0)
        agg = agg.rename(columns={"pIC50_median": "pAct_final", "pIC50_min": "pAct_min",
                                  "pIC50_max": "pAct_max", "pIC50_std": "pAct_std"})
        keep = agg[~agg.discordant]
        first = (grp.sort_values(["pchembl_value", "year"], ascending=False)
                    .drop_duplicates("inchikey").set_index("inchikey"))
        ds = keep.set_index("inchikey").join(
            first[[c for c in PROV if c in first.columns]]).reset_index()
        refs = (grp.groupby("inchikey")
                   .agg(all_assays=("assay_chembl_id", lambda s: "|".join(sorted(set(s.dropna())))),
                        all_docs=("document_chembl_id", lambda s: "|".join(sorted(set(s.dropna())))),
                        all_dois=("doi", lambda s: "|".join(sorted(set(s.dropna())))),
                        all_pmids=("pubmed_id", lambda s: "|".join(
                            sorted(set(str(int(x)) for x in s.dropna()))))))
        ds = ds.merge(refs, on="inchikey", how="left")
        ep = f"{LABEL}_{kind}_{atype}"
        ds["endpoint"] = ep
        ds.to_csv(f"{OUTPUT_DIR}/dataset_{ep}.csv", index=False)
        summary.append({"endpoint": ep, "measure": kind, "assay_type": atype,
                       "records": len(grp), "compounds": grp.inchikey.nunique(),
                       "discarded_discordant": int(agg.discordant.sum()),
                       "compounds_final": len(ds),
                       "pAct_median": round(ds.pAct_final.median(), 2),
                       "frac_actives": round((ds.pAct_final >= 6).mean(), 3),
                       "frac_antagonist": round((ds.assay_mode == "antagonist").mean(), 3),
                       "assays": grp.assay_chembl_id.nunique()})

    summary = pd.DataFrame(summary).sort_values("compounds_final", ascending=False)
    summary.to_csv(f"{OUTPUT_DIR}/summary_pxr.csv", index=False)
    pd.DataFrame(log).to_csv(f"{OUTPUT_DIR}/log_filters_pxr.csv", index=False)
    print(f"\n{len(summary)} endpoints saved to {OUTPUT_DIR}/")
    print(summary.to_string(index=False))
