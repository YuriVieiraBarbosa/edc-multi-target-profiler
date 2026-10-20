# Pipeline QSAR para desreguladores endócrinos

Preditor de consenso multicamada que estima o potencial de moléculas interagirem com receptores
nucleares envolvidos em desregulação endócrina. O projeto cobre a curadoria de dados públicos
(ChEMBL, BindingDB, Tox21), o treino e a validação dos modelos e o uso dos modelos finais para
prever novas moléculas.

## Alvos cobertos

| Receptor | Código |
|---|---|
| Receptor de andrógeno | AR |
| Receptor de estrógeno alfa / beta | ESR1 / ESR2 |
| Receptor do hormônio tireoidiano alfa / beta | THRA / THRB (Tox21: TR, sem distinção alfa/beta) |
| PPAR gama | PPARG |

O receptor de insulina (IR) ainda não é coberto.

## Arquitetura: três camadas independentes

As camadas nunca são fundidas, porque medem coisas diferentes e os ensaios não são comparáveis entre si.

| Camada | Pasta | Modelos | Conteúdo |
|---|---|---|---|
| ChEMBL IC50 | `models_final/` | 18 | regressão de pIC50 e classificação, por alvo e tipo de ensaio (ligação `B` ou funcional `F`) |
| ChEMBL Ki, EC50 e Kd | `models_layers/` | 32 | mesmos alvos, cada medida em seu próprio modelo |
| Tox21 qHTS | `models_tox21/` | 11 | classificação por ensaio (química ambiental), com contra-ensaios de viabilidade |

### Regra de alerta

Para cada par molécula-receptor, só contam as camadas em que a molécula está **dentro do domínio de
aplicabilidade** (similaridade de Tanimoto máxima ao conjunto de treino maior ou igual ao limiar do
modelo) e cujo modelo é implantável. Sobre essas camadas:

| Alerta | Critério |
|---|---|
| alto (farmacologico + ambiental) | 2 ou mais camadas ativas, com pelo menos uma ChEMBL e uma Tox21 |
| moderado (maioria das camadas) | fração de camadas ativas de pelo menos 0,5, com 3 ou mais camadas avaliadas |
| baixo (evidencia isolada) | pelo menos 1 camada ativa |
| sem alerta | nenhuma camada ativa |

Um par em que nenhuma camada está no domínio aparece como `-` na tabela e como `n.d.` nos gráficos.

## Instalação

Requer Python 3.11.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# conferir: deve imprimir 2180
python -c "import features as F; print(2048 + len(F.DESC_NAMES))"
```

As versões de `rdkit`, `scikit-learn` e `xgboost` precisam ser as de `requirements.txt`. Com outra versão
do RDKit o número de descritores muda e o preditor falha com `Feature shape mismatch`.

## Uso

```bash
python rodar_edc.py --input demo_edcs.csv --name-col nome --out demo
```

O CSV de entrada precisa de uma coluna `smiles` (ou informe `--smiles-col`). Cada execução grava em
`resultados/AAAAMMDD_HHMMSS_<nome>/`:

- `<nome>_resumo.csv`: alerta e contagem de camadas por molécula e receptor
- `<nome>_detalhado.csv`: resultado de cada modelo, com domínio, probabilidade ou valor previsto e intervalo conformal
- `<nome>_rejeitadas.csv`: moléculas rejeitadas na padronização, quando houver
- `<nome>_dashboard.png` e `<nome>_fracao_ativa.png`: gráficos de alertas
- `entrada.csv` e `LEIAME.txt`: cópia da entrada, comando, data e versões das bibliotecas

Também é possível usar os scripts separadamente:

```bash
python predict_edc2.py --input demo_edcs.csv --name-col nome --out demo   # só previsão
python plot_edc_results.py --out demo                                     # só gráficos, a partir de demo_*.csv
python predict_edc2.py --smiles "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1"           # uma molécula
```

`painel_ref.csv` é um painel de referência com classes conhecidas e `demo_edcs.csv` é um exemplo curto.

## Conteúdo do repositório

| Caminho | Função |
|---|---|
| `rodar_edc.py` | executa a previsão, gera os gráficos e organiza a saída em `resultados/` |
| `predict_edc2.py` | preditor de consenso multicamada |
| `plot_edc_results.py` | gráficos de alertas e de fração ativa |
| `features.py`, `curation.py` | featurização (ECFP4 por contagem e descritores RDKit) e padronização de estruturas |
| `extract_chembl.py`, `extract_bindingdb.py` | extração de atividades dos dumps do ChEMBL 37 e do BindingDB |
| `pipeline_curate.py`, `pipeline_layers.py`, `pipeline_tox21.py` | curadoria dos datasets de IC50, de Ki/Kd/EC50 e do Tox21 |
| `modeling.py`, `modeling_layers.py`, `modeling_tox21.py` | treino e validação por camada |
| `finalize_models.py`, `conformal.py` | modelos finais e predição conformal com limiar de decisão otimizado |
| `datasets_camadas_por_endpoint/`, `datasets_tox21_por_ensaio/` | datasets curados de Ki/EC50/Kd e de Tox21 usados no treino |
| `models_final/`, `models_layers/`, `models_tox21/` | modelos treinados (cerca de 110 MB) |

Não estão incluídos: os dumps brutos do ChEMBL e do BindingDB, necessários para refazer a extração, e os
datasets curados da camada IC50.

## Protocolo de validação

- Divisão externa por scaffold de Bemis-Murcko (80/20) e seleção de modelo por validação cruzada 5-fold no treino.
- Modelo final reajustado em 100% dos dados curados, com as métricas vindas do teste externo por scaffold.
- Y-scrambling e domínio de aplicabilidade por similaridade de Tanimoto.
- Predição conformal indutiva com cobertura nominal de 80%, e limiar de decisão escolhido por validação cruzada (máximo MCC) em vez de 0,5.
- Curadoria: estruturas válidas, SMILES padronizados, remoção de sais e solventes, detecção de duplicatas e conferência de estereoquímica. Tipos de medida (IC50, Ki, Kd, EC50) e tipos de ensaio nunca são misturados.

## Limitações

- Triagem computacional. Um alerta indica plausibilidade e não substitui ensaio experimental nem avaliação regulatória.
- Moléculas fora do domínio de aplicabilidade não entram na regra de alerta, mesmo quando o modelo devolve uma probabilidade alta. Compostos de classes pouco representadas no treino, como ftalatos, podem ficar sem cobertura.
- Mecanismos que não passam pelos receptores modelados não são detectados. Por exemplo, efeitos via aromatase (CYP19A1) ou esteroidogênese exigiriam endpoints adicionais.
- A camada Tox21 é de química ambiental e usa outro tipo de medida (AC50) que as camadas do ChEMBL.

## Fontes de dados

ChEMBL 37, BindingDB e Tox21. Os dados derivados seguem os termos de cada fonte. O ChEMBL é distribuído sob CC BY-SA 3.0. Consulte os termos das demais fontes antes de redistribuir os datasets.

## Licença

A definir.
