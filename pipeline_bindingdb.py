"""BindingDB layer: the fraction of bioactivity that ChEMBL does not have.

We measured that ChEMBL 37 is exhausted for these targets (loosening confidence_score
adds no records) and that 97.5% of the PubChem BioAssay assays are
deposits from ChEMBL itself. BindingDB is the only substantial increment: 27% of
its records for these targets come from another source, and half of the molecules are
absent from the already curated datasets.

Same funnel and same standardization as the other layers. Two design differences,
both forced by what the source has:

1. PROVENANCE SEPARATES ENDPOINTS. The patent fraction has 86% actives and pAct
   median 7.47, versus 51% and 6.05 from BindingDB's own curation — the same
   enrichment-in-actives signature that led to discarding another source in this
   project. Patent becomes the `_BDBpat` endpoint, separate from `_BDB`, so it can be
   included or excluded downstream without re-curating anything.

2. THERE IS NO ASSAY TYPE OR MODE. BindingDB does not provide an assay description, so
   it's not possible to separate binding from functional as in the ChEMBL layers, nor infer
   agonist/antagonist. The endpoints carry the BDB suffix instead of B/F, and
   `assay_mode` is `undefined` in all rows. It's not an omission: it's absence in the
   source, and needs to travel with any reading of this layer.

Usage:  python pipeline_bindingdb.py [BindingDB_All.tsv] [dataset_IC50_curado.csv]
      The 2nd argument is the ChEMBL IC50 layer, used only to mark which molecules
      are already curated; without it the in_chembl column becomes inflated in the IC50 endpoints.
Output: data_bindingdb/dataset_<ALVO>_<MEDIDA>_BDB[pat].csv, plus the funnel log,
       the per-endpoint summary and the unit check against ChEMBL.
"""
import csv
import glob
import os
import sys

sys.path.insert(0, os.getcwd())
import numpy as np
import pandas as pd

import curation

ACC = {"P10275": "AR", "P03372": "ESR1", "Q92731": "ESR2", "P37231": "PPARG",
       "P10827": "THRA", "P10828": "THRB", "P06213": "IR"}
UP = "UniProt (SwissProt) Primary ID of Target Chain 1"
CHAINS_COLUMN = "Number of Protein Chains in Target (>1 implies a multichain complex)"
ORG = "Target Source Organism According to Curator or DataSource"
MEASURES = ("IC50", "Ki", "Kd", "EC50")
COLS = ["BindingDB Reactant_set_id", "Ligand SMILES", "Ligand InChI Key",
        "BindingDB Ligand Name", "Target Name", ORG,
        "IC50 (nM)", "Ki (nM)", "Kd (nM)", "EC50 (nM)", "pH", "Temp (C)",
        "Curation/DataSource", "Article DOI", "PMID", "PubChem AID", "PubChem CID",
        "ChEMBL ID of Ligand", "Patent Number", "Institution", CHAINS_COLUMN, UP]
MIN_REG = 50
OUTPUT_DIR = "data_bindingdb"

# Columns carried into the final dataset, mirroring the provenance block of the other layers.
PROV = ["bindingdb_id", "ligand_name", "canonical_smiles", "std_smiles", "inchikey",
        "standard_type", "standard_value", "standard_units", "value_nM", "pAct",
        "assay_type", "assay_mode", "assay_organism", "target_name", "accession",
        "target_label", "target_organism", "doi", "pubmed_id", "pubchem_aid",
        "pubchem_cid", "chembl_id_ligand", "patent", "institution", "provenance",
        "source_curation", "source", "undefined_stereo", "mw"]


# Locate the BindingDB TSV: explicit argument, then the usual places in the project tree.
def find_tsv(argument=None):
    """The BindingDB TSV; the 2nd command-line argument is the curated IC50 dataset."""
    candidate_paths = [argument] if argument else []
    candidate_paths += ["bindingdb/BindingDB_All.tsv",
                   os.path.expanduser("~/Downloads/data/bindingdb/BindingDB_All.tsv")]
    for c in candidate_paths:
        if c and os.path.exists(c):
            return c
    raise SystemExit("BindingDB_All.tsv not found. Pass the path as an argument, "
                     "or download it with download_sources.sh.\nSearched in: "
                     + ", ".join(str(c) for c in candidate_paths if c))


# Stream the 8 GB dump and keep only the rows whose target is one of ours.
def extract(tsv):
    chunk_frames, n_read = [], 0
    for ch in pd.read_csv(tsv, sep="\t", usecols=COLS, quoting=csv.QUOTE_NONE,
                          on_bad_lines="skip", low_memory=False, dtype=str,
                          chunksize=200_000):
        n_read += len(ch)
        sel = ch[ch[UP].isin(ACC)]
        if len(sel):
            chunk_frames.append(sel)
    d = pd.concat(chunk_frames, ignore_index=True)
    print(f"rows read from the dump: {n_read} | records in {len(ACC)} targets: {len(d)}", flush=True)
    return d


# One row per (record, measure): a BindingDB row may carry IC50 and Ki at once, and each
# measure becomes its own dataset, exactly as in the ChEMBL layers.
def melt(d):
    output = []
    for m in MEASURES:
        col = f"{m} (nM)"
        sub = d[d[col].notna()].copy()
        if not len(sub):
            continue
        sub["standard_type"] = m
        sub["standard_value"] = sub[col].astype(str).str.strip()
        output.append(sub.drop(columns=[f"{x} (nM)" for x in MEASURES]))
    return pd.concat(output, ignore_index=True)


# Drop censored measurements. The ChEMBL funnel keeps only standard_relation '='; here the
# censoring is inside the value itself (">10000"), and 8-17% of the records carry it.
def without_censoring(d):
    raw_records = d.standard_value.astype(str)
    d = d[~raw_records.str.contains(r"[<>]", na=False)].copy()
    d["value_nM"] = pd.to_numeric(d.standard_value, errors="coerce")
    return d[d.value_nM.notna()]


if __name__ == "__main__":
    tsv = find_tsv(sys.argv[1] if len(sys.argv) > 1 else None)
    print(f"reading: {tsv}", flush=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs("cache", exist_ok=True)

    cache = "cache/bindingdb_targets_raw.pkl"
    raw = pd.read_pickle(cache) if os.path.exists(cache) else extract(tsv)
    if not os.path.exists(cache):
        raw.to_pickle(cache)

    log = []

    # Record the record count after one filtering step, for the curation funnel log.
    def note(step, df):
        log.append({"step": step, "records": len(df)})
        print(f"{step:50s} n={len(df):7d}", flush=True)

    note("0. Raw (4 measures, 7 targets)", raw)
    d = raw[raw[CHAINS_COLUMN] == "1"]
    note("1. single-chain target", d)
    d = d[d[ORG].fillna("").str.contains("Homo sapiens")]
    note("2. human organism", d)
    d = melt(d)
    note("3. one row per record x measure", d)
    d = without_censoring(d)
    note("4. value not censored (without < nor >)", d)
    d = d[d.value_nM.between(1e-3, 1e7)]
    note("5. plausible band (1e-3 to 1e7 nM)", d)
    d = d[d["Ligand SMILES"].notna()]
    note("6. with structure", d)
    # BindingDB ingests ChEMBL: 74% of the records for these targets are reimports of it.
    # Keeping them here would duplicate the datasets that the ChEMBL layers already produce, with the
    # added drawback of passing through another funnel. This layer exists for the COMPLEMENT.
    d = d[d["Curation/DataSource"] != "ChEMBL"]
    note("7. provenance outside ChEMBL", d)

    std_path = "cache/bindingdb_std_lookup.pkl"
    if os.path.exists(std_path):
        std = pd.read_pickle(std_path)
    else:
        std = curation.standardize_frame(d["Ligand SMILES"], n_jobs=10)
        std.to_pickle(std_path)
    d = d.merge(std.rename(columns={"canonical_smiles": "Ligand SMILES"}),
                on="Ligand SMILES", how="left")
    d = d[d.reject_reason.isna()]
    note("8. standardized structure valid", d)

    d["canonical_smiles"] = d["Ligand SMILES"]
    d["standard_units"] = "nM"
    d["pAct"] = curation.to_pic50(d.value_nM)
    d["accession"] = d[UP]
    d["target_label"] = d[UP].map(ACC)
    d["target_organism"] = d[ORG]
    d["assay_organism"] = d[ORG]
    # The source does not provide an assay description: without it there is no way to separate binding from
    # functional nor infer the mode, and inventing a label would be worse than stating the absence.
    d["assay_type"] = "BDB"
    d["assay_mode"] = "undefined"
    d["provenance"] = np.where(d["Curation/DataSource"].eq("US Patent"), "patent", "literature")
    d["source"] = "BindingDB"
    d = d.rename(columns={"BindingDB Reactant_set_id": "bindingdb_id",
                          "BindingDB Ligand Name": "ligand_name",
                          "Target Name": "target_name", "Article DOI": "doi",
                          "PMID": "pubmed_id", "PubChem AID": "pubchem_aid",
                          "PubChem CID": "pubchem_cid",
                          "ChEMBL ID of Ligand": "chembl_id_ligand",
                          "Patent Number": "patent", "Institution": "institution",
                          "Curation/DataSource": "source_curation"})
    d.to_pickle("cache/bindingdb_filtered.pkl")

    # InChIKeys already curated, to mark what is REALLY outside the existing training set.
    # The ChEMBL IC50 datasets may not be in the folder; in that case the mark covers
    # only the present layers, and the script states which ones it found.
    folders = ["data", "ic50_curated_by_endpoint", "datasets_layers_by_endpoint",
              "layers_curated_by_endpoint", "data_layers",
              "datasets_tox21_by_assay", "tox21_curated_by_assay", "data_tox21"]
    found, already = [], set()
    for p in folders:
        fnames = glob.glob(os.path.join(p, "dataset_*.csv"))
        if fnames:
            found.append(f"{p} ({len(fnames)})")
            for a in fnames:
                already |= set(pd.read_csv(a, usecols=["inchikey"]).inchikey.dropna())
    # The ChEMBL IC50 layer is not versioned in the folder. Without it the no_chembl mark
    # ends up INFLATED precisely in the IC50 endpoints, which are the largest in this layer, so
    # its path is passed as the 2nd argument and the absence is stated out loud.
    ic50 = sys.argv[2] if len(sys.argv) > 2 else None
    if ic50 and os.path.exists(ic50):
        already |= set(pd.read_csv(ic50, low_memory=False, usecols=["inchikey"]).inchikey.dropna())
        found.append(f"{os.path.basename(ic50)} (IC50)")
    else:
        print("\nWARNING: without the curated IC50 dataset from ChEMBL (2nd argument). The column "
              "no_chembl is INFLATED in the IC50 endpoints - do not use these numbers this way.")
    print(f"\ndeduplication reference: {len(already)} InChIKeys in {found or 'NO folder'}")
    if not found:
        print("  WARNING: without curated datasets in the folder, the no_chembl column is entirely True")
    d["in_chembl"] = ~d.inchikey.isin(already)

    # ---------------------------------------------------------------- datasets
    summary = []
    for (lab, kind, proc), grp in d.groupby(["target_label", "standard_type", "provenance"]):
        if len(grp) < MIN_REG:
            continue
        agg = curation.aggregate_replicates(grp, ["inchikey"], "pAct", max_spread=1.0)
        agg = agg.rename(columns={"pIC50_median": "pAct_final", "pIC50_min": "pAct_min",
                                  "pIC50_max": "pAct_max", "pIC50_std": "pAct_std"})
        keep = agg[~agg.discordant]
        first = (grp.sort_values("pAct", ascending=False)
                    .drop_duplicates("inchikey").set_index("inchikey"))
        ds = keep.set_index("inchikey").join(
            first[[c for c in PROV if c in first.columns]]).reset_index()
        refs = (grp.groupby("inchikey")
                   .agg(all_dois=("doi", lambda s: "|".join(sorted(set(s.dropna())))),
                        all_pmids=("pubmed_id", lambda s: "|".join(sorted(set(s.dropna())))),
                        all_patents=("patent", lambda s: "|".join(sorted(set(s.dropna())))),
                        sources=("source_curation", lambda s: "|".join(sorted(set(s.dropna()))))))
        ds = ds.merge(refs, on="inchikey", how="left")
        ds["in_chembl"] = ~ds.inchikey.isin(already)
        ep = f"{lab}_{kind}_BDB" + ("pat" if proc == "patent" else "")
        ds["endpoint"] = ep
        ds.to_csv(f"{OUTPUT_DIR}/dataset_{ep}.csv", index=False)
        summary.append({"endpoint": ep, "target": lab, "measure": kind, "provenance": proc,
                       "records": len(grp), "compounds": grp.inchikey.nunique(),
                       "discarded_discordant": int(agg.discordant.sum()),
                       "compounds_final": len(ds),
                       "out_of_chembl": int(ds.in_chembl.sum()),
                       "frac_out_of_chembl": round(ds.in_chembl.mean(), 3),
                       "pAct_median": round(ds.pAct_final.median(), 2),
                       "frac_actives": round((ds.pAct_final >= 6).mean(), 3)})

    summary = pd.DataFrame(summary).sort_values("compounds_final", ascending=False)
    summary.to_csv(f"{OUTPUT_DIR}/summary_bindingdb.csv", index=False)
    pd.DataFrame(log).to_csv(f"{OUTPUT_DIR}/log_filters_bindingdb.csv", index=False)
    print(f"\n{len(summary)} endpoints saved to {OUTPUT_DIR}/")
    print(summary.to_string(index=False))
