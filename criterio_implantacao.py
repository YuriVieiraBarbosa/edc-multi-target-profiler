#!/usr/bin/env python
"""Criterio explicito de implantacao dos modelos.

O flag `nao_implantavel` gravado dentro dos .joblib de models_final/ veio de um script que
nao esta no repositorio, entao nao era auditavel: `ESR2_B` passava com R2 de 0,098 e `AR_F`
era reprovado com 0,076, sem regra escrita em lugar nenhum. Este script reconstroi o
criterio a partir das metricas que os proprios artefatos guardam, grava
`criterio_implantacao.csv` e passa a ser a fonte que o predict_edc2.py consulta.

Regra, aplicada por MODELO (nao por endpoint):

    n_train < 200                      -> reprovado, em qualquer camada
    regressao:    Q2_cv < 0,30         -> reprovado
                  R2_teste_scaffold<0  -> reprovado (pior que prever a media)
    classificacao: ROC_AUC_teste < 0,60 -> reprovado

Tox21 e um caso a parte: nem os .joblib nem o LEIAME guardam metrica de validacao, entao
esses 11 modelos NAO sao avaliaveis por este criterio. Eles seguem votando, mas por
ausencia de criterio, nao por terem passado em um - o que esta registrado na coluna
`motivo` e deve ser dito em qualquer relatorio.

Uso: python criterio_implantacao.py [saida.csv]
"""
import glob
import os
import sys

import joblib
import numpy as np
import pandas as pd

N_MINIMO = 200
Q2_MINIMO = 0.30
R2_TESTE_MINIMO = 0.0
ROC_AUC_TESTE_MINIMO = 0.60


PASTAS = ("models_final", "models_layers", "models_tox21")
ARQUIVO = "models_reprovados"     # reprovados guardados so como documentacao (ver arquivar_reprovados.py)


# Yield (logical folder, path, where) for every trained .joblib: the run folders and the archive.
# The logical folder stays the ORIGINAL one, because the rule and the filename grammar depend on it.
def _modelos(raiz="."):
    for pasta in PASTAS:
        for local, d in (("execucao", os.path.join(raiz, pasta)),
                         ("arquivo", os.path.join(raiz, ARQUIVO, pasta))):
            for p in sorted(glob.glob(os.path.join(d, "*.joblib"))):
                yield pasta, p, local


# Read every trained model (running or archived) and collect the metrics each layer stores.
def coletar(raiz="."):
    cam = os.path.join(raiz, "models_layers", "LEIAME_modelos_camadas.csv")
    lay = pd.read_csv(cam).set_index("endpoint") if os.path.exists(cam) else pd.DataFrame()
    linhas = []
    for pasta, p, local in _modelos(raiz):
        m = joblib.load(p)
        base = os.path.basename(p)[:-7]
        # o campo 'task' nao e consistente entre camadas (ha 'regression' e 'regressao'),
        # entao a tarefa vem do sufixo do nome do arquivo, que e uniforme nos modelos
        reg = base.endswith("_regressao")
        r = dict(pasta=pasta, modelo_id=base, local=local,
                 tarefa="regressao" if reg else "classificacao",
                 flag_gravado=bool(m.get("nao_implantavel", False)),
                 n=np.nan, metrica_cv=np.nan, metrica_teste=np.nan)
        if pasta == "models_final":
            r["n"] = m.get("n_train")
            r["metrica_cv"] = m.get("Q2_cv") if reg else m.get("ROC_AUC_cv")
            r["metrica_teste"] = (m.get("R2_teste_scaffold") if reg
                                  else m.get("ROC_AUC_teste_scaffold"))
        elif pasta == "models_layers":
            ep = base.rsplit("_", 1)[0]
            if ep in lay.index:
                r["n"] = lay.loc[ep, "n"]
                if reg:
                    r["metrica_cv"] = lay.loc[ep, "Q2_cv"]
        else:
            r["n"] = m.get("n_train")
        linhas.append(r)
    return pd.DataFrame(linhas)


# Apply the rule above to one model, returning (rejected?, reason).
def julgar(r):
    if r.pasta == "models_tox21":
        return False, "sem metrica de validacao gravada (nao avaliavel)"
    if pd.notna(r.n) and r.n < N_MINIMO:
        return True, f"n={int(r.n)} < {N_MINIMO}"
    if r.tarefa == "regressao":
        if pd.notna(r.metrica_cv) and r.metrica_cv < Q2_MINIMO:
            return True, f"Q2_cv={r.metrica_cv:.3f} < {Q2_MINIMO:.2f}"
        if pd.notna(r.metrica_teste) and r.metrica_teste < R2_TESTE_MINIMO:
            return True, f"R2_teste_scaffold={r.metrica_teste:.3f} < {R2_TESTE_MINIMO:.2f}"
    elif pd.notna(r.metrica_teste) and r.metrica_teste < ROC_AUC_TESTE_MINIMO:
        return True, f"ROC_AUC_teste={r.metrica_teste:.3f} < {ROC_AUC_TESTE_MINIMO:.2f}"
    return False, "aprovado"


# Build the full table: one row per model, with the verdict and the reason behind it.
def tabela(raiz="."):
    T = coletar(raiz)
    T[["nao_implantavel", "motivo"]] = T.apply(lambda r: pd.Series(julgar(r)), axis=1)
    T["divergiu_do_flag_gravado"] = T.nao_implantavel != T.flag_gravado
    # a ordem nao pode depender de onde o arquivo esta hoje, senao mover um modelo
    # reescreve o CSV inteiro e o diff do git deixa de mostrar o que mudou de fato
    T = T.sort_values(["pasta", "modelo_id"]).reset_index(drop=True)
    cols = [c for c in T.columns if c != "local"] + ["local"]
    return T[cols]


if __name__ == "__main__":
    saida = sys.argv[1] if len(sys.argv) > 1 else "criterio_implantacao.csv"
    T = tabela()
    T.to_csv(saida, index=False)
    n_rep = int(T.nao_implantavel.sum())
    print(f"{len(T)} modelos | implantados: {len(T) - n_rep} | reprovados: {n_rep}")
    desloc = T[(T.local == "execucao") & T.nao_implantavel]
    if len(desloc):
        print(f"  ATENCAO: {len(desloc)} reprovado(s) ainda nas pastas de execucao; "
              "rode arquivar_reprovados.py --aplicar")
    print(f"  dos implantados, {int((T.pasta == 'models_tox21').sum())} sao do Tox21 e nao "
          "foram avaliados por falta de metrica gravada")
    div = T[T.divergiu_do_flag_gravado]
    if len(div):
        print(f"\ndivergencias do flag gravado no .joblib ({len(div)}):")
        for _, r in div.iterrows():
            print(f"  {r.modelo_id:28s} gravado={'reprovado' if r.flag_gravado else 'aprovado':9s}"
                  f" criterio={'reprovado' if r.nao_implantavel else 'aprovado':9s} ({r.motivo})")
    print(f"\ngravado em {saida}")
