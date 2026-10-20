"""Camadas paralelas de bioatividade do ChEMBL 37: Ki, Kd e EC50.

Mesmos alvos, mesmos filtros de comparabilidade e mesma padronização do
dataset de IC50. Cada tipo de medida gera datasets independentes —
Ki, Kd, EC50 e IC50 nunca são fundidos.
"""
import os, re, sys, sqlite3
sys.path.append(os.getcwd())
import numpy as np
import pandas as pd
import curation
from extract_chembl import ACT_SQL

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
        "doi", "pubmed_id", "year", "journal", "activity_source", "source",
        "undefined_stereo", "mw"]

DB = "chembl/chembl_37/chembl_37_sqlite/chembl_37.db"
TIPOS = ["Ki", "Kd", "EC50"]
MIN_REG = 50

if not os.path.exists("cache/chembl37_targets.csv"):
    from extract_chembl import TARGET_SQL, ACC
    _con = sqlite3.connect(DB)
    _tg = pd.read_sql(TARGET_SQL.format(placeholders=",".join("?" * len(ACC))),
                      _con, params=list(ACC.values()))
    _con.close()
    _tg["target_label"] = _tg["accession"].map({v: k for k, v in ACC.items()})
    os.makedirs("cache", exist_ok=True)
    _tg.to_csv("cache/chembl37_targets.csv", index=False)
tg = pd.read_csv("cache/chembl37_targets.csv")
single = tg[(tg.target_type == "SINGLE PROTEIN") & (tg.organism == "Homo sapiens")]
tids = sorted(set(single.tid))

con = sqlite3.connect(DB)
sql = ACT_SQL.format(tids=",".join("?" * len(tids))).replace(
    "AND act.standard_type = 'IC50'",
    "AND act.standard_type IN ({})".format(",".join("?" * len(TIPOS))))
raw = pd.read_sql(sql, con, params=tids + TIPOS)
con.close()
raw = raw.merge(single[["target_chembl_id", "accession", "target_label"]].drop_duplicates(),
                on="target_chembl_id", how="left")

log = []
# Record the record count after one filtering step, for the curation funnel log.
def note(step, df):
    log.append({"etapa": step, "registros": len(df)})
    print(f"{step:50s} n={len(df):7d}")

note("0. brutos (Ki/Kd/EC50, alvo humano)", raw)
d = raw[raw.target_organism == "Homo sapiens"]
d = d[d.confidence_score >= 8]; note("1. confidence_score >= 8", d)
d = d[d.standard_relation == "="]; note("2. relação '='", d)
d = d[d.standard_units.isin(UNIT_TO_NM)].copy()
d["value_nM"] = pd.to_numeric(d.standard_value, errors="coerce") * d.standard_units.map(UNIT_TO_NM)
note("3. unidade conversível para nM", d)
d = d[d.value_nM.between(1e-3, 1e7)]; note("4. faixa plausível", d)
d = d[~d.data_validity_comment.isin(BAD_VALIDITY)]; note("5. validade ok", d)
d = d[d.canonical_smiles.notna() & (d.structure_type == "MOL")]; note("6. com estrutura", d)

std_path = "cache/layers_std_lookup.pkl"
if os.path.exists(std_path):
    std = pd.read_pickle(std_path)
else:
    std = curation.standardize_frame(d.canonical_smiles, n_jobs=10)
    std.to_pickle(std_path)
d = d.merge(std, on="canonical_smiles", how="left")
d = d[d.reject_reason.isna()]; note("7. estrutura padronizada válida", d)

d["pAct"] = curation.to_pic50(d.value_nM)
d["assay_mode"] = np.where(d.assay_description.fillna("").str.contains(ANTAG_RE), "antagonista",
                  np.where(d.assay_description.fillna("").str.contains(AGONIST_RE), "agonista", "indefinido"))
d["source"] = "ChEMBL37"
d.to_pickle("cache/chembl_layers_filtrado.pkl")

os.makedirs("data_layers", exist_ok=True)
resumo = []
for (lab, tipo, atype), grp in d.groupby(["target_label", "standard_type", "assay_type"]):
    if atype not in ("B", "F") or len(grp) < MIN_REG:
        continue
    agg = curation.aggregate_replicates(grp, ["inchikey"], "pAct", max_spread=1.0)
    agg = agg.rename(columns={"pIC50_median": "pAct_final", "pIC50_min": "pAct_min",
                              "pIC50_max": "pAct_max", "pIC50_std": "pAct_std"})
    keep = agg[~agg.discordant]
    first = (grp.sort_values(["pchembl_value", "year"], ascending=False)
                .drop_duplicates("inchikey").set_index("inchikey"))
    ds = keep.set_index("inchikey").join(first[[c for c in PROV if c in first.columns]]).reset_index()
    refs = (grp.groupby("inchikey")
               .agg(all_assays=("assay_chembl_id", lambda s: "|".join(sorted(set(s.dropna())))),
                    all_docs=("document_chembl_id", lambda s: "|".join(sorted(set(s.dropna())))),
                    all_dois=("doi", lambda s: "|".join(sorted(set(s.dropna())))),
                    all_pmids=("pubmed_id", lambda s: "|".join(sorted(set(str(int(x)) for x in s.dropna()))))))
    ds = ds.merge(refs, on="inchikey", how="left")
    ep = f"{lab}_{tipo}_{atype}"
    ds["endpoint"] = ep
    ds.to_csv(f"data_layers/dataset_{ep}.csv", index=False)
    resumo.append({"endpoint": ep, "alvo": lab, "medida": tipo, "assay_type": atype,
                   "registros": len(grp), "compostos": grp.inchikey.nunique(),
                   "descartados_discordantes": int(agg.discordant.sum()),
                   "compostos_finais": len(ds),
                   "pAct_mediana": round(ds.pAct_final.median(), 2),
                   "frac_antagonista": round((ds.assay_mode == "antagonista").mean(), 3)})

resumo = pd.DataFrame(resumo).sort_values("compostos_finais", ascending=False)
resumo.to_csv("data_layers/resumo_camadas.csv", index=False)
pd.DataFrame(log).to_csv("data_layers/log_filtros_camadas.csv", index=False)
print("\n", resumo.to_string(index=False))
