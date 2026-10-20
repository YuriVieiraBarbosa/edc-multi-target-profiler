#!/usr/bin/env python
"""Consenso v2: relatório de dois eixos, sem retreinar nada.

Lê o `<out>_detalhado.csv` produzido por predict_edc2.py e reescreve só a agregação. Nenhum modelo é
alterado ou retreinado, e os votos individuais são exatamente os mesmos.

Mudanças em relação ao `_resumo.csv` atual:
  1. A unidade de evidência é o ENDPOINT, não o modelo. A regressão e a classificação do mesmo endpoint
     ChEMBL deixam de contar duas vezes: a classificação decide a positividade (generaliza melhor fora
     do treino) e a regressão entra como ordem de grandeza da potência, com o intervalo conformal.
  2. Dois eixos, NUNCA somados: farmacológico (ChEMBL) e triagem (Tox21). Cada um tem o seu estado.
  3. Dentro do eixo, a MELHOR evidência e quantos endpoints a corroboram, em vez da fração de votos
     (os endpoints não são independentes: Ki e EC50 compartilham compostos, e os ensaios Tox21 de um
     receptor são os mesmos modos em linhagens diferentes).
  4. Três estados por eixo: positivo, negativo e `nao_avaliado`. A cobertura é explícita
     (endpoints que existem / que puderam votar), e quem está fora do domínio nunca dilui nada.
  5. A direção do efeito viaja como rótulo (do modo dos ensaios Tox21), fora do escore.
  6. Toda combinação molécula x receptor ganha linha, inclusive sem nenhuma cobertura. No formato atual
     esses pares simplesmente não apareciam, e a ausência de linha era confundida com "sem alerta".

Uso:
    python consenso_v2.py <out>_detalhado.csv [prefixo_de_saida]

Saídas:
    <prefixo>_endpoints.csv  uma linha por molécula x receptor x endpoint (unidade de evidência)
    <prefixo>_resumo.csv     uma linha por molécula x receptor, com os dois eixos e a cobertura
    <prefixo>_cobertura.csv  endpoints disponíveis por receptor e eixo (o "teto" de cada receptor)
"""
import sys

import numpy as np
import pandas as pd

EIXO = {"ChEMBL": "farmacologico", "Tox21": "triagem"}
BOOL = ("conta", "ativo", "no_dominio", "nao_implantavel")


# Read the predictor's _detalhado.csv and normalise its boolean, optional and derived columns.
def carregar(path):
    d = pd.read_csv(path)
    for c in BOOL:
        d[c] = d[c].astype(str).str.lower().eq("true")
    # colunas opcionais: o predict_edc2.py só cria intervalo_inf/sup quando há calibração conformal,
    # e prob_ativo/limiar só quando há classificador. Sem isto, um arquivo sem conformal quebra abaixo.
    for c in ("valor_previsto", "intervalo_inf", "intervalo_sup", "prob_ativo", "limiar", "tanimoto_viz"):
        if c not in d.columns:
            d[c] = np.nan
    if "detalhe" not in d.columns:
        d["detalhe"] = ""
    d["eixo"] = np.where(d.camada == "Tox21", "triagem", "farmacologico")
    # endpoint base: o nome sem o sufixo da tarefa (só as camadas ChEMBL têm as duas tarefas)
    d["endpoint_base"] = d.endpoint.str.replace(r"_(regressao|classificacao)$", "", regex=True)
    d["medida"] = np.where(d.eixo == "triagem", "AC50", d.camada.str.replace("ChEMBL ", "", regex=False))
    return d


# Build one evidence row per voting unit, for every molecule x receptor x axis.
def por_endpoint(d, unidade="modelo"):
    """Uma linha de evidência por unidade de voto.

    unidade="modelo"   (padrão): cada modelo vota sozinho. Regressão e classificação do
                       mesmo endpoint são DOIS votos independentes. É o esquema original
                       do predict_edc2.py.
    unidade="endpoint": colapsa regressão + classificação do mesmo endpoint em UMA
                       evidência (a classificação decide, a regressão vira potência).
    """
    if unidade not in ("modelo", "endpoint"):
        raise ValueError(f"unidade deve ser 'modelo' ou 'endpoint', recebido {unidade!r}")
    chave = ["nome", "receptor", "eixo", "medida", "endpoint_base"]
    if unidade == "modelo":
        # incluir o endpoint com sufixo de tarefa deixa cada modelo num grupo só seu
        chave = chave + ["endpoint"]
    linhas = []
    for k, g in d.groupby(chave, sort=True):      # sort=True: a ordem das linhas não depende da entrada
        cls = g[g.tarefa == "classificacao"]
        reg = g[g.tarefa == "regressao"]
        for nome_t, t in (("classificacao", cls), ("regressao", reg)):
            if len(t) > 1:                            # nunca escolher em silêncio qual linha vale
                raise ValueError(
                    f"{len(t)} linhas de {nome_t} para o mesmo endpoint em {dict(zip(chave, k))}. "
                    "Causa habitual: duas moléculas com o mesmo valor na coluna de nome (ou SMILES "
                    "repetido sem --name-col), que o predict_edc2.py funde numa linha só.")
        impl = bool((~g.nao_implantavel).any())          # existe algum modelo implantável neste endpoint
        cls_c, reg_c = cls[cls.conta], reg[reg.conta]
        if len(cls_c):                                    # a classificação decide a positividade
            base, positivo = "classificacao", bool(cls_c.ativo.iloc[0])
            prob, limiar = float(cls_c.prob_ativo.iloc[0]), float(cls_c.limiar.iloc[0])
        elif len(reg_c):                                  # sem classificador no domínio: cai para a regressão
            base, positivo = "regressao", bool(reg_c.ativo.iloc[0])
            prob, limiar = np.nan, np.nan
        else:
            base, positivo, prob, limiar = "nao_avaliado", None, np.nan, np.nan
        r = dict(zip(chave, k))
        r.update(implantavel=impl, avaliado=base != "nao_avaliado", base_decisao=base,
                 positivo=positivo, prob_ativo=prob, limiar=limiar,
                 tanimoto_viz=float(g.tanimoto_viz.max()),
                 modo=g.detalhe.iloc[0] if r["eixo"] == "triagem" else "",
                 pAct_previsto=float(reg_c.valor_previsto.iloc[0]) if len(reg_c) else np.nan,
                 pAct_inf=float(reg_c.intervalo_inf.iloc[0]) if len(reg_c) else np.nan,
                 pAct_sup=float(reg_c.intervalo_sup.iloc[0]) if len(reg_c) else np.nan)
        linhas.append(r)
    e = pd.DataFrame(linhas)
    e["modo_simples"] = e.modo.str.split("/").str[0]
    return e


# State of one axis for a molecule x receptor pair, with vote counts, coverage and best evidence.
def _estado(g):
    """Estado de um eixo para um par molécula x receptor, a partir dos endpoints daquele eixo."""
    disp = g[g.implantavel]
    av = disp[disp.avaliado]
    pos = av[av.positivo == True]
    if len(disp) == 0:
        est = "sem_modelo"
    elif len(av) == 0:
        est = "nao_avaliado"
    else:
        est = "positivo" if len(pos) else "negativo"
    # melhor evidência: maior probabilidade entre os endpoints positivos (ou entre os avaliados)
    ref = pos if len(pos) else av
    i = ref.prob_ativo.idxmax() if len(ref) and ref.prob_ativo.notna().any() else (ref.index[0] if len(ref) else None)
    melhor = ref.loc[i] if i is not None else None
    # Potência: pIC50, pKi e pEC50 NÃO são a mesma escala, então o máximo entre medidas não é
    # interpretável sozinho. Guardamos de qual endpoint e medida ele veio, e também a potência de cada
    # medida separadamente, que é a leitura correta quando o receptor tem mais de uma.
    pot = av[av.pAct_previsto.notna()]
    j = pot.pAct_previsto.idxmax() if len(pot) else None
    por_medida = "; ".join(f"{med} {pot.loc[idx, 'pAct_previsto']:.2f}"
                           for med, idx in sorted(pot.groupby("medida").pAct_previsto.idxmax().items()))
    # "votos recebidos / quantidade de modelos", nas duas leituras possíveis:
    #   fracao_positiva           -> denominador = unidades que votaram (esquema anterior)
    #   fracao_sobre_implantaveis -> denominador = todas as unidades implantáveis do eixo
    # A diferença entre as duas é a cobertura do domínio de aplicabilidade.
    return pd.Series(dict(
        estado=est, endpoints_disponiveis=len(disp), endpoints_avaliados=len(av), endpoints_positivos=len(pos),
        fracao_positiva=(round(len(pos) / len(av), 2) if len(av) else np.nan),
        fracao_sobre_implantaveis=(round(len(pos) / len(disp), 2) if len(disp) else np.nan),
        melhor_endpoint=(melhor.endpoint_base if melhor is not None else ""),
        melhor_prob=(melhor.prob_ativo if melhor is not None else np.nan),
        melhor_base=(melhor.base_decisao if melhor is not None else ""),
        pAct_melhor=(pot.loc[j, "pAct_previsto"] if j is not None else np.nan),
        pAct_melhor_inf=(pot.loc[j, "pAct_inf"] if j is not None else np.nan),
        pAct_melhor_sup=(pot.loc[j, "pAct_sup"] if j is not None else np.nan),
        pAct_melhor_endpoint=(pot.loc[j, "endpoint_base"] if j is not None else ""),
        pAct_melhor_medida=(pot.loc[j, "medida"] if j is not None else ""),
        potencia_por_medida=por_medida,
        tanimoto_max=(disp.tanimoto_viz.max() if len(disp) else np.nan),
        medidas_positivas=", ".join(sorted(pos.medida.unique())),
        modos_positivos=", ".join(sorted(set(pos.modo_simples.dropna()) - {""})),
    ))


# One row per molecule x receptor, with both axes side by side and the joint reading.
def resumir(e):
    # indexa uma vez por (nome, receptor, eixo) em vez de refiltrar o frame a cada par
    blocos = {k: g for k, g in e.groupby(["nome", "receptor", "eixo"], sort=False)}
    vazio = e.iloc[0:0]
    pares = e[["nome", "receptor"]].drop_duplicates().sort_values(["nome", "receptor"])
    out = []
    for _, (nome, rec) in pares.iterrows():
        r = dict(nome=nome, receptor=rec)
        for eixo, pref in (("farmacologico", "farm"), ("triagem", "tox")):
            s = _estado(blocos.get((nome, rec, eixo), vazio))
            r.update({f"{pref}_{k}": v for k, v in s.items()})
        out.append(r)
    r = pd.DataFrame(out)
    # leitura conjunta dos dois eixos, sem somar nada: todas as 9 combinações de estado são nomeadas
    LEITURA = {
        ("positivo", "positivo"): "convergente positivo",
        ("positivo", "negativo"): "divergente: farmacologico positivo, triagem negativo",
        ("positivo", "sem cobertura"): "só farmacologico positivo (sem cobertura de triagem)",
        ("negativo", "positivo"): "divergente: triagem positivo, farmacologico negativo",
        ("sem cobertura", "positivo"): "só triagem positivo (sem cobertura farmacologica)",
        ("negativo", "negativo"): "convergente negativo",
        ("negativo", "sem cobertura"): "negativo só no eixo farmacologico (cobertura parcial)",
        ("sem cobertura", "negativo"): "negativo só no eixo de triagem (cobertura parcial)",
        ("sem cobertura", "sem cobertura"): "sem cobertura",
    }
    simp = lambda s: np.where(s.isin(["nao_avaliado", "sem_modelo"]), "sem cobertura", s)
    r["leitura"] = [LEITURA[(a, b)] for a, b in zip(simp(r.farm_estado), simp(r.tox_estado))]
    # direção: só o Tox21 distingue modo; o ChEMBL não carrega agonista/antagonista
    r["direcao"] = np.where(r.tox_modos_positivos == "", "nao determinada", r.tox_modos_positivos)
    cols = (["nome", "receptor", "leitura", "direcao"]
            + [c for c in r.columns if c.startswith("farm_")] + [c for c in r.columns if c.startswith("tox_")])
    return r[cols]


# Voting ceiling per receptor and axis: how many units could vote at most.
def cobertura(e):
    """Unidades implantáveis por receptor e eixo: o teto de votos, igual para toda molécula."""
    # a unidade é o modelo quando por_endpoint rodou com unidade="modelo"
    unid = "endpoint" if "endpoint" in e.columns else "endpoint_base"
    d = e[e.implantavel].drop_duplicates(["receptor", "eixo", unid])
    t = d.pivot_table(index="receptor", columns="eixo", values=unid, aggfunc="count", fill_value=0)
    for c in ("farmacologico", "triagem"):
        if c not in t:
            t[c] = 0
    t["medidas_farmacologicas"] = (d[d.eixo == "farmacologico"].groupby("receptor").medida
                                   .apply(lambda s: ", ".join(sorted(set(s)))))
    t["modos_triagem"] = (d[d.eixo == "triagem"].groupby("receptor").modo_simples
                          .apply(lambda s: ", ".join(sorted(set(s)))))
    return t.fillna("").reset_index()


if __name__ == "__main__":
    # uso: consenso_v2.py <detalhado.csv> [prefixo_saida] [modelo|endpoint]
    det = sys.argv[1]
    pref = sys.argv[2] if len(sys.argv) > 2 else "consenso_v2"
    unidade = sys.argv[3] if len(sys.argv) > 3 else "modelo"
    d = carregar(det)
    e = por_endpoint(d, unidade=unidade)
    r = resumir(e)
    c = cobertura(e)
    e.to_csv(f"{pref}_endpoints.csv", index=False)
    r.to_csv(f"{pref}_resumo.csv", index=False)
    c.to_csv(f"{pref}_cobertura.csv", index=False)
    print(f"unidade de voto: {unidade}")
    print(f"{len(d)} previsões -> {len(e)} unidades -> {len(r)} pares molécula x receptor")
    print(r.leitura.value_counts().to_string())
