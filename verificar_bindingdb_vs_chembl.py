"""Confere a camada BindingDB contra o ChEMBL nos compostos que as duas fontes têm.

Uma camada nova só é confiável depois de bater com a fonte já curada nos pontos em
comum. Neste projeto essa conferência já pagou: numa fonte anterior ela revelou um
grupo de pares que divergiam por um fator de mil, concentrado em poucas origens, com
assinatura de erro de unidade.

O que procura, por par (composto, alvo, medida) presente nas duas fontes:
  |dif| <= 1,0    dentro da tolerância de réplicas que o próprio funil usa
  |dif| > 1,0     divergência real, listada
  |dif| ~ 3,0     candidato a erro de unidade (fator 1.000, tipicamente nM x uM)

Uso:  python verificar_bindingdb_vs_chembl.py [caminho/do/dataset_IC50_curado.csv]
Saída: data_bindingdb/conferencia_vs_chembl.csv e o resumo na tela.

A camada de IC50 do ChEMBL não está versionada na pasta; passe o caminho dela como
argumento para incluí-la, senão a conferência cobre só Ki, Kd e EC50.
"""
import glob
import os
import sys

sys.path.insert(0, os.getcwd())
import numpy as np
import pandas as pd
from scipy import stats

SAIDA = "data_bindingdb"


# Load the curated ChEMBL side: the Ki/Kd/EC50 layers plus, if given, the IC50 dataset.
def lado_chembl(ic50_path=None):
    partes = []
    pastas = ["datasets_camadas_por_endpoint", "camadas_curado_por_endpoint", "data_layers"]
    for p in pastas:
        for a in glob.glob(os.path.join(p, "dataset_*.csv")):
            d = pd.read_csv(a, usecols=["inchikey", "endpoint", "pAct_final"])
            d["alvo"] = d.endpoint.str.split("_").str[0]
            d["medida"] = d.endpoint.str.split("_").str[1]
            partes.append(d[["inchikey", "alvo", "medida", "pAct_final"]])
        if partes:
            print(f"camadas ChEMBL lidas de: {p}")
            break
    if ic50_path and os.path.exists(ic50_path):
        d = pd.read_csv(ic50_path, low_memory=False, usecols=["inchikey", "endpoint", "pIC50_final"])
        d["alvo"] = d.endpoint.str.split("_").str[0]
        d["medida"] = "IC50"
        partes.append(d.rename(columns={"pIC50_final": "pAct_final"})[
            ["inchikey", "alvo", "medida", "pAct_final"]])
        print(f"camada de IC50 lida de: {ic50_path}")
    else:
        print("AVISO: sem o dataset de IC50 do ChEMBL; a conferencia cobre so Ki, Kd e EC50")
    if not partes:
        raise SystemExit("nenhum dataset curado do ChEMBL encontrado")
    return pd.concat(partes, ignore_index=True)


# Load the curated BindingDB side written by pipeline_bindingdb.py.
def lado_bindingdb():
    arqs = sorted(glob.glob(f"{SAIDA}/dataset_*_BDB*.csv"))
    if not arqs:
        raise SystemExit(f"nenhum dataset em {SAIDA}/. Rode pipeline_bindingdb.py antes.")
    partes = []
    for a in arqs:
        d = pd.read_csv(a, low_memory=False,
                        usecols=["inchikey", "endpoint", "pAct_final", "target_label",
                                 "standard_type", "procedencia", "fontes"])
        partes.append(d)
    return pd.concat(partes, ignore_index=True)


if __name__ == "__main__":
    ch = lado_chembl(sys.argv[1] if len(sys.argv) > 1 else None)
    bd = lado_bindingdb()
    print(f"\npares curados: ChEMBL {len(ch)} | BindingDB {len(bd)}")

    m = bd.merge(ch, left_on=["inchikey", "target_label", "standard_type"],
                 right_on=["inchikey", "alvo", "medida"],
                 suffixes=("_bdb", "_chembl"))
    m = (m.groupby(["inchikey", "target_label", "standard_type", "procedencia", "fontes"],
                   as_index=False)
           .agg(pAct_bdb=("pAct_final_bdb", "median"),
                pAct_chembl=("pAct_final_chembl", "median")))
    m["dif"] = m.pAct_bdb - m.pAct_chembl
    m["suspeita_unidade"] = m.dif.abs().between(2.5, 3.5)
    m.to_csv(f"{SAIDA}/conferencia_vs_chembl.csv", index=False)

    if not len(m):
        print("\nnenhum par em comum: nada a conferir")
        sys.exit(0)

    print(f"\n=== {len(m)} pares (composto, alvo, medida) nas duas fontes ===")
    print(f"  diferenca mediana BindingDB - ChEMBL : {m.dif.median():+.2f} log")
    print(f"  |dif| mediana                        : {m.dif.abs().median():.2f} log")
    print(f"  dentro da tolerancia de 1,0 log      : {int((m.dif.abs() <= 1).sum())} "
          f"({(m.dif.abs() <= 1).mean():.1%})")
    print(f"  |dif| > 1,0                          : {int((m.dif.abs() > 1).sum())}")
    print(f"  |dif| entre 2,5 e 3,5 (fator ~1.000) : {int(m.suspeita_unidade.sum())}")
    if len(m) > 2:
        print(f"  Spearman entre as duas fontes         : {stats.spearmanr(m.pAct_bdb, m.pAct_chembl).statistic:.2f}")

    print("\n=== por procedencia ===")
    print(m.groupby("procedencia").agg(
        pares=("dif", "size"), dif_mediana=("dif", "median"),
        absdif_mediana=("dif", lambda s: s.abs().median()),
        acima_de_1=("dif", lambda s: int((s.abs() > 1).sum())),
        suspeita_unidade=("suspeita_unidade", "sum")).round(2).to_string())

    sus = m[m.suspeita_unidade]
    if len(sus):
        print(f"\n=== {len(sus)} pares com assinatura de erro de unidade ===")
        print(sus.nlargest(min(12, len(sus)), "dif", keep="all")[
            ["inchikey", "target_label", "standard_type", "procedencia",
             "pAct_bdb", "pAct_chembl", "dif", "fontes"]].round(2).to_string(index=False))
    print(f"\ngravado: {SAIDA}/conferencia_vs_chembl.csv")
