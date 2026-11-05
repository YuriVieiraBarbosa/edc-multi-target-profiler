#!/usr/bin/env python
"""Mantem os modelos reprovados fora da execucao, guardados so como documentacao.

Um modelo reprovado por criterio_implantacao.csv nunca vota. Deixa-lo em models_final/,
models_layers/ ou models_tox21/ so custa tempo: o preditor o carrega e o aplica a toda
molecula para depois descartar o resultado. Este script move os reprovados para
models_reprovados/<pasta de origem>/, de onde o preditor nao carrega nada, e mantem o
registro de cada um em models_reprovados/LEIAME_reprovados.csv.

E simetrico: um modelo arquivado que passe a ser APROVADO (por exemplo depois de um
retreino com mais dados e de rodar criterio_implantacao.py de novo) volta para a pasta de
execucao. Modelo nenhum e apagado ou alterado, so trocado de pasta.

Fluxo depois de treinar modelos novos:
    python criterio_implantacao.py          # regrava o criterio, vendo as duas pastas
    python arquivar_reprovados.py           # mostra o que mudaria de pasta
    python arquivar_reprovados.py --aplicar # move de fato

Uso: python arquivar_reprovados.py [--aplicar]
"""
import os
import shutil
import sys

import pandas as pd

ARQUIVO = "models_reprovados"
CRITERIO = "criterio_implantacao.csv"

LEIAME = """# Modelos reprovados (so documentacao)

Estes modelos foram treinados e **nao sao carregados pelo preditor**: `predict_edc2.py`, `edc.py`
e `rodar_edc.py` olham apenas `models_final/`, `models_layers/` e `models_tox21/`. Estao aqui para
registrar o que foi tentado e por que nao entrou. Cada um tem a linha correspondente em
`LEIAME_reprovados.csv` e em `../criterio_implantacao.csv`.

"Reprovado" quase sempre quer dizer **pouco dado**, nao modelo ruim: o criterio exige n >= 200 e a
maior parte dos reprovados falha so nisso. O corte de 200 e uma convencao do projeto, nao um valor
validado. Veja a coluna `motivo`.

Os indices `models_final/LEIAME_calibracao.csv` e `models_layers/LEIAME_modelos_camadas.csv`
continuam listando todos os endpoints treinados, inclusive estes, de modo que um endpoint listado
la pode nao ter `.joblib` na pasta de execucao.

## Como um modelo volta para a execucao

Nao ha mover a mao. Retreine com mais dados, rode `python criterio_implantacao.py` e depois
`python arquivar_reprovados.py --aplicar`: o que passar no criterio volta para a pasta de origem e o
que reprovar vem para ca. O arquivo `.joblib` mantem o mesmo nome nos dois lugares.
"""


# Compare where each model is with where the criterion says it belongs; return the moves.
def planejar(raiz="."):
    t = pd.read_csv(os.path.join(raiz, CRITERIO))
    mover = []
    for r in t.itertuples():
        exe = os.path.join(raiz, r.pasta, f"{r.modelo_id}.joblib")
        arq = os.path.join(raiz, ARQUIVO, r.pasta, f"{r.modelo_id}.joblib")
        if r.nao_implantavel and os.path.exists(exe):
            mover.append((r.modelo_id, "arquivar", exe, arq, r.motivo))
        elif not r.nao_implantavel and os.path.exists(arq):
            mover.append((r.modelo_id, "promover", arq, exe, r.motivo))
    return t, mover


# Write the documentation of the archived models, from the criterion table.
def documentar(t, raiz="."):
    rep = t[t.nao_implantavel][["modelo_id", "pasta", "tarefa", "n", "metrica_cv",
                                "metrica_teste", "motivo"]]
    rep = rep.rename(columns={"pasta": "pasta_de_origem"})
    os.makedirs(os.path.join(raiz, ARQUIVO), exist_ok=True)
    rep.to_csv(os.path.join(raiz, ARQUIVO, "LEIAME_reprovados.csv"), index=False)
    with open(os.path.join(raiz, ARQUIVO, "LEIAME.md"), "w") as fh:
        fh.write(LEIAME)
    return len(rep)


if __name__ == "__main__":
    aplicar = "--aplicar" in sys.argv[1:]
    if not os.path.exists(CRITERIO):
        raise SystemExit(f"{CRITERIO} nao encontrado. Rode python criterio_implantacao.py antes.")
    t, mover = planejar()
    n_arq = sum(1 for m in mover if m[1] == "arquivar")
    n_pro = len(mover) - n_arq
    print(f"{len(t)} modelos no criterio | a arquivar: {n_arq} | a promover: {n_pro}")
    conflito = [m for m in mover if os.path.exists(m[3])]
    if conflito:
        for m in conflito:
            print(f"  CONFLITO: {m[0]} existe nos dois lugares; nada foi movido.")
        raise SystemExit(1)
    for nome, acao, origem, destino, motivo in mover:
        print(f"  {acao:8s} {nome:28s} ({motivo})")
    if not aplicar:
        print("\nsimulacao: nada foi movido. Use --aplicar para mover.")
    else:
        for nome, acao, origem, destino, motivo in mover:
            os.makedirs(os.path.dirname(destino), exist_ok=True)
            shutil.move(origem, destino)
        n = documentar(t)
        print(f"\nmovidos {len(mover)} | {n} modelos documentados em {ARQUIVO}/LEIAME_reprovados.csv")
