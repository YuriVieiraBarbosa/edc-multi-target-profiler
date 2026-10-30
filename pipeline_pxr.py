"""Camada ChEMBL do receptor de pregnano X (PXR / NR1I2, UniProt O75469).

O PXR entra no projeto como xenossensor: ele regula a CYP3A4 e medeia boa parte do
efeito dos desreguladores por inducao enzimatica, alem de ter ligantes em comum com
ER, AR e PPARgama. Mesmo funil, mesma padronizacao e mesma gramatica de endpoint das
outras camadas do ChEMBL.

UMA DIFERENCA DE DESENHO, forcada pelo que a fonte tem. No PXR o campo `assay_type` do
ChEMBL separa mal ligacao de funcional - defeito ja registrado no projeto para os
outros alvos, mas aqui ele e extremo:

  tipo A (ADMET)  421 registros, o MAIOR grupo: na descricao sao reporter de CYP3A4 em
                  HepG2/DPX-2, ou seja FUNCIONAIS, mais um TR-FRET de ligacao
  tipo B          646 registros: mistura deslocamento de BODIPY FL vindoline (ligacao
                  de verdade) com reporter de luciferase em HepG2 (funcional)
  tipo F           66 registros, reporter
  tipo T           32 registros, transativacao

Descartar o tipo A, como o funil das outras camadas faz, jogaria fora 36% dos registros
do alvo. Fundi-lo ao B violaria a regra do projeto de nunca misturar ensaios
incompativeis. A saida adotada e a terceira: cada assay_type vira seu PROPRIO endpoint
(PXR_EC50_A, PXR_IC50_B, ...), que e a mesma regra ja aplicada a B e F, so que estendida.
A consequencia tem de viajar com qualquer leitura: nenhum endpoint de PXR pode ser
chamado de 'ligacao' ou 'funcional' pelo sufixo.

Uso:  python pipeline_pxr.py [caminho/do/chembl_37.db]
Saida: data_pxr/dataset_PXR_<MEDIDA>_<TIPO>.csv, log do funil e resumo por endpoint.
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
ROTULO = "PXR"
TIPOS = ["IC50", "Ki", "Kd", "EC50"]
MIN_REG = 50
SAIDA = "data_pxr"

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
def achar_db(argumento=None):
    cand = [argumento] if argumento else []
    cand += ["chembl/chembl_37/chembl_37_sqlite/chembl_37.db",
             os.path.expanduser("~/Downloads/dados/chembl_37/chembl_37/"
                                "chembl_37_sqlite/chembl_37.db")]
    for c in cand:
        if c and os.path.exists(c):
            return c
    raise SystemExit("chembl_37.db nao encontrado. Passe o caminho como argumento.\n"
                     "Procurei em: " + ", ".join(str(c) for c in cand if c))


if __name__ == "__main__":
    db = achar_db(sys.argv[1] if len(sys.argv) > 1 else None)
    print(f"lendo: {db}", flush=True)
    os.makedirs(SAIDA, exist_ok=True)
    os.makedirs("cache", exist_ok=True)

    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    tg = pd.read_sql(TARGET_SQL.format(placeholders="?"), con, params=[ACC_PXR])
    single = tg[(tg.target_type == "SINGLE PROTEIN") & (tg.organism == "Homo sapiens")]
    if not len(single):
        raise SystemExit(f"nenhum alvo SINGLE PROTEIN humano para {ACC_PXR}")
    tids = sorted(set(single.tid))
    print(f"alvo: {single.target_chembl_id.iloc[0]} "
          f"({single.target_name.iloc[0]}) | tid {tids}", flush=True)

    sql = ACT_SQL.format(tids=",".join("?" * len(tids))).replace(
        "AND act.standard_type = 'IC50'",
        "AND act.standard_type IN ({})".format(",".join("?" * len(TIPOS))))
    raw = pd.read_sql(sql, con, params=tids + TIPOS)
    con.close()
    raw = raw.merge(single[["target_chembl_id", "accession"]].drop_duplicates(),
                    on="target_chembl_id", how="left")
    raw["target_label"] = ROTULO

    log = []

    # Record the record count after one filtering step, for the curation funnel log.
    def note(etapa, df):
        log.append({"etapa": etapa, "registros": len(df)})
        print(f"{etapa:50s} n={len(df):7d}", flush=True)

    note("0. brutos (IC50/Ki/Kd/EC50, PXR humano)", raw)
    d = raw[raw.target_organism == "Homo sapiens"]
    d = d[d.confidence_score >= 8]
    note("1. confidence_score >= 8", d)
    d = d[d.standard_relation == "="]
    note("2. relacao '='", d)
    d = d[d.standard_units.isin(UNIT_TO_NM)].copy()
    d["value_nM"] = pd.to_numeric(d.standard_value, errors="coerce") * d.standard_units.map(UNIT_TO_NM)
    note("3. unidade conversivel para nM", d)
    d = d[d.value_nM.between(1e-3, 1e7)]
    note("4. faixa plausivel", d)
    d = d[~d.data_validity_comment.isin(BAD_VALIDITY)]
    note("5. validade ok", d)
    d = d[d.canonical_smiles.notna() & (d.structure_type == "MOL")]
    note("6. com estrutura", d)

    std_path = "cache/pxr_std_lookup.pkl"
    if os.path.exists(std_path):
        std = pd.read_pickle(std_path)
    else:
        std = curation.standardize_frame(d.canonical_smiles, n_jobs=10)
        std.to_pickle(std_path)
    d = d.merge(std, on="canonical_smiles", how="left")
    d = d[d.reject_reason.isna()]
    note("7. estrutura padronizada valida", d)

    d["pAct"] = curation.to_pic50(d.value_nM)
    # O modo vem da descricao, como nas outras camadas. A ordem da regex testa antagonista
    # primeiro, entao 'inhibit'/'block' ganham de 'agonis' - limitacao ja conhecida.
    desc = d.assay_description.fillna("")
    d["assay_mode"] = np.where(desc.str.contains(ANTAG_RE), "antagonista",
                      np.where(desc.str.contains(AGONIST_RE), "agonista", "indefinido"))
    d["source"] = "ChEMBL37"
    d.to_pickle("cache/pxr_filtrado.pkl")

    resumo = []
    for (tipo, atype), grp in d.groupby(["standard_type", "assay_type"]):
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
        ep = f"{ROTULO}_{tipo}_{atype}"
        ds["endpoint"] = ep
        ds.to_csv(f"{SAIDA}/dataset_{ep}.csv", index=False)
        resumo.append({"endpoint": ep, "medida": tipo, "assay_type": atype,
                       "registros": len(grp), "compostos": grp.inchikey.nunique(),
                       "descartados_discordantes": int(agg.discordant.sum()),
                       "compostos_finais": len(ds),
                       "pAct_mediana": round(ds.pAct_final.median(), 2),
                       "frac_ativos": round((ds.pAct_final >= 6).mean(), 3),
                       "frac_antagonista": round((ds.assay_mode == "antagonista").mean(), 3),
                       "ensaios": grp.assay_chembl_id.nunique()})

    resumo = pd.DataFrame(resumo).sort_values("compostos_finais", ascending=False)
    resumo.to_csv(f"{SAIDA}/resumo_pxr.csv", index=False)
    pd.DataFrame(log).to_csv(f"{SAIDA}/log_filtros_pxr.csv", index=False)
    print(f"\n{len(resumo)} endpoints gravados em {SAIDA}/")
    print(resumo.to_string(index=False))
