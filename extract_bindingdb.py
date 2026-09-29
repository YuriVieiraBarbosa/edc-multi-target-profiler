"""Filtra o dump completo do BindingDB (TSV) pelas UniProt dos alvos de
desregulação endócrina, mantendo apenas registros com IC50."""
import csv, pandas as pd

ACC = {"P10275": "AR", "P03372": "ESR1", "Q92731": "ESR2",
       "P37231": "PPARG", "P10827": "THRA", "P10828": "THRB"}

COLS = ["BindingDB Reactant_set_id", "Ligand SMILES", "Ligand InChI Key",
        "BindingDB Ligand Name", "Target Name",
        "Target Source Organism According to Curator or DataSource",
        "IC50 (nM)", "Ki (nM)", "Kd (nM)", "EC50 (nM)", "pH", "Temp (C)",
        "Curation/DataSource", "Article DOI", "PMID", "PubChem AID", "PubChem CID",
        "ChEMBL ID of Ligand", "Patent Number", "Institution",
        "Number of Protein Chains in Target (>1 implies a multichain complex)",
        "UniProt (SwissProt) Primary ID of Target Chain 1",
        "UniProt (SwissProt) Entry Name of Target Chain 1"]

UP = "UniProt (SwissProt) Primary ID of Target Chain 1"

out, n, kept = [], 0, 0
reader = pd.read_csv("bindingdb/BindingDB_All.tsv", sep="\t", usecols=COLS,
                     quoting=csv.QUOTE_NONE, on_bad_lines="skip", low_memory=False,
                     dtype=str, chunksize=200_000)
for ch in reader:
    n += len(ch)
    sel = ch[ch[UP].isin(ACC) & ch["IC50 (nM)"].notna()]
    if len(sel):
        out.append(sel)
        kept += len(sel)

df = pd.concat(out, ignore_index=True)
df["target_label"] = df[UP].map(ACC)
df.to_parquet("cache/bindingdb_ic50_raw.parquet", index=False)
print("linhas lidas:", n, "| registros IC50 nos alvos:", kept)
print(df.groupby("target_label").size().to_string())
print(df["Number of Protein Chains in Target (>1 implies a multichain complex)"].value_counts().to_dict())
