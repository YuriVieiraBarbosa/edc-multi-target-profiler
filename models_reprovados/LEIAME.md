# Modelos reprovados (so documentacao)

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
