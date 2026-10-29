"""Camada BindingDB: a fração de bioatividade que o ChEMBL não tem.

Medimos que o ChEMBL 37 está esgotado para estes alvos (afrouxar o confidence_score
não acrescenta nenhum registro) e que 97,5% dos ensaios do PubChem BioAssay são
depósito do próprio ChEMBL. O BindingDB é o único incremento substancial: 27% dos
seus registros nestes alvos vêm de outra procedência, e metade das moléculas está
ausente dos datasets já curados.

Mesmo funil e mesma padronização das outras camadas. Duas diferenças de desenho,
ambas forçadas pelo que a fonte tem:

1. PROCEDÊNCIA SEPARA ENDPOINTS. A fração de patentes tem 86% de ativos e pAct
   mediana 7,47, contra 51% e 6,05 da curadoria do próprio BindingDB — a mesma
   assinatura de enriquecimento em ativos que levou a descartar outra fonte neste
   projeto. Patente vira endpoint `_BDBpat`, separado do `_BDB`, para poder ser
   incluída ou excluída a jusante sem recurar nada.

2. NÃO HÁ TIPO DE ENSAIO NEM MODO. O BindingDB não traz descrição de ensaio, então
   não dá para separar ligação de funcional como nas camadas do ChEMBL, nem inferir
   agonista/antagonista. Os endpoints levam o sufixo BDB em vez de B/F, e
   `assay_mode` fica `indefinido` em todas as linhas. Não é omissão: é ausência na
   fonte, e precisa viajar com qualquer leitura desta camada.

Uso:  python pipeline_bindingdb.py [BindingDB_All.tsv] [dataset_IC50_curado.csv]
      O 2o argumento e a camada de IC50 do ChEMBL, usada so para marcar quais moleculas
      ja estao curadas; sem ele a coluna no_chembl fica inflada nos endpoints de IC50.
Saída: data_bindingdb/dataset_<ALVO>_<MEDIDA>_BDB[pat].csv, mais o log do funil,
       o resumo por endpoint e a conferência de unidade contra o ChEMBL.
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
CADEIAS = "Number of Protein Chains in Target (>1 implies a multichain complex)"
ORG = "Target Source Organism According to Curator or DataSource"
MEDIDAS = ("IC50", "Ki", "Kd", "EC50")
COLS = ["BindingDB Reactant_set_id", "Ligand SMILES", "Ligand InChI Key",
        "BindingDB Ligand Name", "Target Name", ORG,
        "IC50 (nM)", "Ki (nM)", "Kd (nM)", "EC50 (nM)", "pH", "Temp (C)",
        "Curation/DataSource", "Article DOI", "PMID", "PubChem AID", "PubChem CID",
        "ChEMBL ID of Ligand", "Patent Number", "Institution", CADEIAS, UP]
MIN_REG = 50
SAIDA = "data_bindingdb"

# Columns carried into the final dataset, mirroring the provenance block of the other layers.
PROV = ["bindingdb_id", "ligand_name", "canonical_smiles", "std_smiles", "inchikey",
        "standard_type", "standard_value", "standard_units", "value_nM", "pAct",
        "assay_type", "assay_mode", "assay_organism", "target_name", "accession",
        "target_label", "target_organism", "doi", "pubmed_id", "pubchem_aid",
        "pubchem_cid", "chembl_id_ligante", "patente", "instituicao", "procedencia",
        "fonte_curadoria", "source", "undefined_stereo", "mw"]


# Locate the BindingDB TSV: explicit argument, then the usual places in the project tree.
def achar_tsv(argumento=None):
    """O TSV do BindingDB; o 2o argumento da linha de comando e o dataset de IC50 curado."""
    candidatos = [argumento] if argumento else []
    candidatos += ["bindingdb/BindingDB_All.tsv",
                   os.path.expanduser("~/Downloads/dados/bindingdb/BindingDB_All.tsv")]
    for c in candidatos:
        if c and os.path.exists(c):
            return c
    raise SystemExit("BindingDB_All.tsv nao encontrado. Passe o caminho como argumento, "
                     "ou baixe com baixar_fontes.sh.\nProcurei em: "
                     + ", ".join(str(c) for c in candidatos if c))


# Stream the 8 GB dump and keep only the rows whose target is one of ours.
def extrair(tsv):
    pedacos, lidas = [], 0
    for ch in pd.read_csv(tsv, sep="\t", usecols=COLS, quoting=csv.QUOTE_NONE,
                          on_bad_lines="skip", low_memory=False, dtype=str,
                          chunksize=200_000):
        lidas += len(ch)
        sel = ch[ch[UP].isin(ACC)]
        if len(sel):
            pedacos.append(sel)
    d = pd.concat(pedacos, ignore_index=True)
    print(f"linhas lidas do dump: {lidas} | registros nos {len(ACC)} alvos: {len(d)}", flush=True)
    return d


# One row per (record, measure): a BindingDB row may carry IC50 and Ki at once, and each
# measure becomes its own dataset, exactly as in the ChEMBL layers.
def derreter(d):
    saida = []
    for m in MEDIDAS:
        col = f"{m} (nM)"
        sub = d[d[col].notna()].copy()
        if not len(sub):
            continue
        sub["standard_type"] = m
        sub["standard_value"] = sub[col].astype(str).str.strip()
        saida.append(sub.drop(columns=[f"{x} (nM)" for x in MEDIDAS]))
    return pd.concat(saida, ignore_index=True)


# Drop censored measurements. The ChEMBL funnel keeps only standard_relation '='; here the
# censoring is inside the value itself (">10000"), and 8-17% of the records carry it.
def sem_censura(d):
    bruto = d.standard_value.astype(str)
    d = d[~bruto.str.contains(r"[<>]", na=False)].copy()
    d["value_nM"] = pd.to_numeric(d.standard_value, errors="coerce")
    return d[d.value_nM.notna()]


if __name__ == "__main__":
    tsv = achar_tsv(sys.argv[1] if len(sys.argv) > 1 else None)
    print(f"lendo: {tsv}", flush=True)
    os.makedirs(SAIDA, exist_ok=True)
    os.makedirs("cache", exist_ok=True)

    cache = "cache/bindingdb_alvos_raw.pkl"
    raw = pd.read_pickle(cache) if os.path.exists(cache) else extrair(tsv)
    if not os.path.exists(cache):
        raw.to_pickle(cache)

    log = []

    # Record the record count after one filtering step, for the curation funnel log.
    def note(etapa, df):
        log.append({"etapa": etapa, "registros": len(df)})
        print(f"{etapa:50s} n={len(df):7d}", flush=True)

    note("0. brutos (4 medidas, 7 alvos)", raw)
    d = raw[raw[CADEIAS] == "1"]
    note("1. alvo de cadeia unica", d)
    d = d[d[ORG].fillna("").str.contains("Homo sapiens")]
    note("2. organismo humano", d)
    d = derreter(d)
    note("3. uma linha por registro x medida", d)
    d = sem_censura(d)
    note("4. valor nao censurado (sem < nem >)", d)
    d = d[d.value_nM.between(1e-3, 1e7)]
    note("5. faixa plausivel (1e-3 a 1e7 nM)", d)
    d = d[d["Ligand SMILES"].notna()]
    note("6. com estrutura", d)
    # O BindingDB ingere o ChEMBL: 74% dos registros destes alvos sao reimportacao dele.
    # Mante-los aqui duplicaria os datasets que as camadas do ChEMBL ja produzem, com a
    # agravante de passarem por outro funil. Esta camada existe para o COMPLEMENTO.
    d = d[d["Curation/DataSource"] != "ChEMBL"]
    note("7. procedencia fora do ChEMBL", d)

    std_path = "cache/bindingdb_std_lookup.pkl"
    if os.path.exists(std_path):
        std = pd.read_pickle(std_path)
    else:
        std = curation.standardize_frame(d["Ligand SMILES"], n_jobs=10)
        std.to_pickle(std_path)
    d = d.merge(std.rename(columns={"canonical_smiles": "Ligand SMILES"}),
                on="Ligand SMILES", how="left")
    d = d[d.reject_reason.isna()]
    note("8. estrutura padronizada valida", d)

    d["canonical_smiles"] = d["Ligand SMILES"]
    d["standard_units"] = "nM"
    d["pAct"] = curation.to_pic50(d.value_nM)
    d["accession"] = d[UP]
    d["target_label"] = d[UP].map(ACC)
    d["target_organism"] = d[ORG]
    d["assay_organism"] = d[ORG]
    # A fonte nao traz descricao de ensaio: sem ela nao ha como separar ligacao de
    # funcional nem inferir o modo, e inventar um rotulo seria pior que declarar a ausencia.
    d["assay_type"] = "BDB"
    d["assay_mode"] = "indefinido"
    d["procedencia"] = np.where(d["Curation/DataSource"].eq("US Patent"), "patente", "literatura")
    d["source"] = "BindingDB"
    d = d.rename(columns={"BindingDB Reactant_set_id": "bindingdb_id",
                          "BindingDB Ligand Name": "ligand_name",
                          "Target Name": "target_name", "Article DOI": "doi",
                          "PMID": "pubmed_id", "PubChem AID": "pubchem_aid",
                          "PubChem CID": "pubchem_cid",
                          "ChEMBL ID of Ligand": "chembl_id_ligante",
                          "Patent Number": "patente", "Institution": "instituicao",
                          "Curation/DataSource": "fonte_curadoria"})
    d.to_pickle("cache/bindingdb_filtrado.pkl")

    # InChIKeys ja curados, para marcar o que esta REALMENTE fora do treino existente.
    # Os datasets de IC50 do ChEMBL podem nao estar na pasta; nesse caso a marca cobre
    # so as camadas presentes, e o script diz quais encontrou.
    pastas = ["data", "ic50_curado_por_endpoint", "datasets_camadas_por_endpoint",
              "camadas_curado_por_endpoint", "data_layers",
              "datasets_tox21_por_ensaio", "tox21_curado_por_ensaio", "data_tox21"]
    achadas, ja = [], set()
    for p in pastas:
        arqs = glob.glob(os.path.join(p, "dataset_*.csv"))
        if arqs:
            achadas.append(f"{p} ({len(arqs)})")
            for a in arqs:
                ja |= set(pd.read_csv(a, usecols=["inchikey"]).inchikey.dropna())
    # A camada de IC50 do ChEMBL nao esta versionada na pasta. Sem ela a marca no_chembl
    # fica INFLADA justamente nos endpoints de IC50, que sao os maiores desta camada, entao
    # o caminho dela entra como 2o argumento e a ausencia e dita em voz alta.
    ic50 = sys.argv[2] if len(sys.argv) > 2 else None
    if ic50 and os.path.exists(ic50):
        ja |= set(pd.read_csv(ic50, low_memory=False, usecols=["inchikey"]).inchikey.dropna())
        achadas.append(f"{os.path.basename(ic50)} (IC50)")
    else:
        print("\nAVISO: sem o dataset de IC50 curado do ChEMBL (2o argumento). A coluna "
              "no_chembl fica INFLADA nos endpoints de IC50 - nao use esses numeros assim.")
    print(f"\nreferencia de deduplicacao: {len(ja)} InChIKeys em {achadas or 'NENHUMA pasta'}")
    if not achadas:
        print("  AVISO: sem datasets curados na pasta, a coluna no_chembl fica toda True")
    d["no_chembl"] = ~d.inchikey.isin(ja)

    # ---------------------------------------------------------------- datasets
    resumo = []
    for (lab, tipo, proc), grp in d.groupby(["target_label", "standard_type", "procedencia"]):
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
                        all_patentes=("patente", lambda s: "|".join(sorted(set(s.dropna())))),
                        fontes=("fonte_curadoria", lambda s: "|".join(sorted(set(s.dropna()))))))
        ds = ds.merge(refs, on="inchikey", how="left")
        ds["no_chembl"] = ~ds.inchikey.isin(ja)
        ep = f"{lab}_{tipo}_BDB" + ("pat" if proc == "patente" else "")
        ds["endpoint"] = ep
        ds.to_csv(f"{SAIDA}/dataset_{ep}.csv", index=False)
        resumo.append({"endpoint": ep, "alvo": lab, "medida": tipo, "procedencia": proc,
                       "registros": len(grp), "compostos": grp.inchikey.nunique(),
                       "descartados_discordantes": int(agg.discordant.sum()),
                       "compostos_finais": len(ds),
                       "fora_do_chembl": int(ds.no_chembl.sum()),
                       "frac_fora_do_chembl": round(ds.no_chembl.mean(), 3),
                       "pAct_mediana": round(ds.pAct_final.median(), 2),
                       "frac_ativos": round((ds.pAct_final >= 6).mean(), 3)})

    resumo = pd.DataFrame(resumo).sort_values("compostos_finais", ascending=False)
    resumo.to_csv(f"{SAIDA}/resumo_bindingdb.csv", index=False)
    pd.DataFrame(log).to_csv(f"{SAIDA}/log_filtros_bindingdb.csv", index=False)
    print(f"\n{len(resumo)} endpoints gravados em {SAIDA}/")
    print(resumo.to_string(index=False))
