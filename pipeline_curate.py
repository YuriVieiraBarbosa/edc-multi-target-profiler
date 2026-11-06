"""Curation pipeline: ChEMBL 37 (primary) + BindingDB (secondary/external).

Produces, per endpoint and per assay kind, comparable pIC50 datasets,
with full provenance and a log of counts per filter.
"""
import os, sys, re, json
sys.path.append(os.getcwd())
import numpy as np
import pandas as pd
import curation

os.makedirs("data", exist_ok=True)

BAD_VALIDITY = {"Outside typical range", "Non standard unit for type",
                "Potential missing data", "Potential transcription error",
                "Potential author error", "Author confirmed error",
                "Potential duplicate"}
UNIT_TO_NM = {"nM": 1.0, "uM": 1e3, "µM": 1e3, "mM": 1e6, "M": 1e9, "pM": 1e-3}

AGONIST_RE = re.compile(r"\bagonis", re.I)
ANTAG_RE = re.compile(r"\bantagonis|\binhibit|\bblock", re.I)

log = []


# Record the record and compound count after one filtering step, for the curation funnel log.
def note(step, df, extra=""):
    log.append({"step": step, "records": len(df),
                "compounds": df["molecule_chembl_id"].nunique() if "molecule_chembl_id" in df else np.nan,
                "obs": extra})
    print(f"{step:55s} n={len(df):7d} {extra}")


# ---------------------------------------------------------------- ChEMBL
act = pd.read_parquet("cache/chembl37_ic50_raw.parquet")
note("0. Raw IC50 (human target SINGLE PROTEIN)", act)

act = act[act["target_organism"] == "Homo sapiens"]
note("1. Homo sapiens target", act)

act = act[act["confidence_score"] >= 8]
note("2. confidence_score >= 8 (direct protein target)", act)

act = act[act["standard_relation"] == "="]
note("3. standard_relation '='", act)

act = act[act["standard_units"].isin(UNIT_TO_NM)]
act["value_nM"] = pd.to_numeric(act["standard_value"], errors="coerce") * \
    act["standard_units"].map(UNIT_TO_NM)
note("4. unit convertible to nM", act)

act = act[act["value_nM"].between(1e-3, 1e7, inclusive="both")]
note("5. value in plausible band (pIC50 2-12)", act)

act = act[~act["data_validity_comment"].isin(BAD_VALIDITY)]
note("6. without problematic data_validity_comment", act)

act = act[act["canonical_smiles"].notna() & (act["structure_type"] == "MOL")]
note("7. with defined chemical structure", act)

# structural standardization
std_path = "cache/chembl_std_lookup.pkl"
if os.path.exists(std_path):
    std = pd.read_pickle(std_path)
else:
    std = curation.standardize_frame(act["canonical_smiles"], n_jobs=10)
    std.to_pickle(std_path)
print("  SMILES unique standardized:", len(std),
      "| rejected:", int(std["reject_reason"].notna().sum()))
rej = std[std["reject_reason"].notna()]["reject_reason"].str.split(":").str[0].value_counts()
print("  reasons:", rej.to_dict())

act = act.merge(std, on="canonical_smiles", how="left")
act_rejected = act[act["reject_reason"].notna()].copy()
act = act[act["reject_reason"].isna()]
note("8. valid standardized structure (salt/solvent removed)", act)

act["pIC50"] = curation.to_pic50(act["value_nM"])
act["assay_mode"] = np.where(
    act["assay_description"].fillna("").str.contains(ANTAG_RE), "antagonist",
    np.where(act["assay_description"].fillna("").str.contains(AGONIST_RE), "agonist", "undefined"))
act["source"] = "ChEMBL37"
act.to_pickle("cache/chembl_ic50_filtered.pkl")

# ---------------------------------------------------------------- BindingDB
# extract_bindingdb.py writes .parquet; .pkl remains accepted for old caches.
bdb_raw = "cache/bindingdb_ic50_raw.parquet"
if os.path.exists(bdb_raw):
    bdb = pd.read_parquet(bdb_raw)
elif os.path.exists("cache/bindingdb_ic50_raw.pkl"):
    bdb = pd.read_pickle("cache/bindingdb_ic50_raw.pkl")
else:
    raise SystemExit(f"{bdb_raw} not found. Run before: python extract_bindingdb.py")
bdb = bdb[bdb["Number of Protein Chains in Target (>1 implies a multichain complex)"] == "1"]
raw_ic = bdb["IC50 (nM)"].astype(str).str.strip()
bdb = bdb[~raw_ic.str.contains(r"[<>]")]
bdb["value_nM"] = pd.to_numeric(raw_ic.str.replace(r"[^0-9eE.+-]", "", regex=True),
                                errors="coerce")
bdb = bdb[bdb["value_nM"].between(1e-3, 1e7)]
# No script wrote cache/bindingdb_std_lookup.pkl, which before was only READ here.
# The lookup now gets built like ChEMBL's (same function, same cache on disk),
# and after the numeric filters, to standardize fewer SMILES.
bstd_path = "cache/bindingdb_std_lookup.pkl"
if os.path.exists(bstd_path):
    bstd = pd.read_pickle(bstd_path)
else:
    bstd = curation.standardize_frame(bdb["Ligand SMILES"], n_jobs=10)
    bstd.to_pickle(bstd_path)
print("  BindingDB SMILES unique standardized:", len(bstd),
      "| rejected:", int(bstd["reject_reason"].notna().sum()))
bstd = bstd.rename(columns={"canonical_smiles": "Ligand SMILES"})
bdb = bdb.merge(bstd, on="Ligand SMILES", how="left")
bdb = bdb[bdb["reject_reason"].isna()]
bdb["pIC50"] = curation.to_pic50(bdb["value_nM"])
bdb["source"] = "BindingDB_202609"
print(f"BindingDB after equivalent filters: n={len(bdb)}")
bdb.to_pickle("cache/bindingdb_ic50_filtered.pkl")

# ------------------------------------------------- aggregation by endpoint
PROV = ["molecule_chembl_id", "molecule_name", "canonical_smiles", "std_smiles",
        "inchikey", "standard_type", "standard_relation", "standard_value",
        "standard_units", "value_nM", "pIC50", "pchembl_value", "assay_chembl_id",
        "assay_description", "assay_type", "assay_mode", "assay_organism",
        "confidence_score", "bao_format", "target_chembl_id", "target_name",
        "accession", "target_label", "target_organism", "document_chembl_id",
        "doi", "pubmed_id", "year", "journal", "activity_source", "source",
        "undefined_stereo", "mw", "activity_id"]

records, datasets, summary = [], {}, []
for (lab, atype), grp in act.groupby(["target_label", "assay_type"]):
    if atype not in ("B", "F") or len(grp) < 50:
        continue
    agg = curation.aggregate_replicates(grp, ["inchikey"], "pIC50", max_spread=1.0)
    keep = agg[~agg["discordant"]]
    first = (grp.sort_values(["pchembl_value", "year"], ascending=False)
                .drop_duplicates("inchikey").set_index("inchikey"))
    ds = keep.set_index("inchikey").join(first[[c for c in PROV if c in first.columns]])
    ds = ds.reset_index().rename(columns={"pIC50_median": "pIC50_final"})
    ds["endpoint"] = f"{lab}_{atype}"
    # provenance of all references that support each compound
    refs = (grp.groupby("inchikey")
               .agg(all_assays=("assay_chembl_id", lambda s: "|".join(sorted(set(s.dropna())))),
                    all_docs=("document_chembl_id", lambda s: "|".join(sorted(set(s.dropna())))),
                    all_dois=("doi", lambda s: "|".join(sorted(set(s.dropna())))),
                    all_pmids=("pubmed_id", lambda s: "|".join(sorted(set(str(int(x)) for x in s.dropna()))))))
    ds = ds.merge(refs, on="inchikey", how="left")
    datasets[f"{lab}_{atype}"] = ds
    summary.append({"endpoint": f"{lab}_{atype}", "target": lab, "assay_type": atype,
                    "records": len(grp), "compounds_unique": grp["inchikey"].nunique(),
                    "discarded_discordant": int(agg["discordant"].sum()),
                    "compounds_final": len(ds),
                    "pIC50_min": round(ds["pIC50_final"].min(), 2),
                    "pIC50_max": round(ds["pIC50_final"].max(), 2),
                    "pIC50_median": round(ds["pIC50_final"].median(), 2),
                    "frac_antagonist_mode": round((ds["assay_mode"] == "antagonist").mean(), 3)})
    ds.to_csv(f"data/dataset_{lab}_{atype}.csv", index=False)

summary = pd.DataFrame(summary).sort_values("compounds_final", ascending=False)
summary.to_csv("data/endpoints_summary.csv", index=False)
pd.DataFrame(log).to_csv("data/log_filters_chembl.csv", index=False)
print("\n", summary.to_string(index=False))

# ------------------------------------- independent external set (BindingDB)
ext_rows = []
for lab, grp in bdb.groupby("target_label"):
    chembl_keys = set(act[act["target_label"] == lab]["inchikey"])
    g = grp[~grp["inchikey"].isin(chembl_keys)]
    agg = curation.aggregate_replicates(g, ["inchikey"], "pIC50", max_spread=1.0)
    keep = agg[~agg["discordant"]].rename(columns={"pIC50_median": "pIC50_final"})
    first = g.drop_duplicates("inchikey").set_index("inchikey")
    e = keep.set_index("inchikey").join(first[["std_smiles", "Ligand SMILES", "BindingDB Ligand Name",
                                               "Target Name", "Article DOI", "PMID",
                                               "Curation/DataSource", "ChEMBL ID of Ligand",
                                               "PubChem CID", "source"]]).reset_index()
    e["target_label"] = lab
    ext_rows.append(e)
ext = pd.concat(ext_rows, ignore_index=True)
ext.to_csv("data/external_bindingdb.csv", index=False)
print("\nExternal BindingDB set (compounds absent from ChEMBL):")
print(ext.groupby("target_label").size().to_string())
