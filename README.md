# EDCQSAR Endocrine Disruption Chemical Quantitative Structure-Activity Relationship

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
| Receptor de pregnano X | PXR (NR1I2) |

O receptor de insulina (IR) ainda não é coberto.

O **PXR** entrou depois dos demais, por `pipeline_pxr.py` (ChEMBL) e `pipeline_pxr_tox21.py` (Tox21).
Ele não é um receptor hormonal como os outros: é o xenossensor que regula a CYP3A4, e entra porque
medeia parte do efeito dos desreguladores por indução enzimática e compartilha ligantes com ER, AR e
PPARγ. Três particularidades dele, todas medidas:

- **O sufixo de tipo de ensaio não indica o que o ensaio mede.** O campo `assay_type` do ChEMBL separa
  mal ligação de funcional em todos os alvos, mas no PXR isso é extremo: o maior grupo funcional
  (repórter de CYP3A4 em HepG2/DPX-2) vem rotulado `A` (ADMET), e o `B` mistura deslocamento de BODIPY
  FL vindoline, que é ligação de verdade, com repórter de luciferase. Descartar o `A`, como o funil das
  outras camadas faz, jogaria fora 36% dos registros do alvo; fundi-lo ao `B` violaria a regra de nunca
  misturar ensaios incompatíveis. A saída adotada foi a terceira: **cada `assay_type` é seu próprio
  endpoint**, que é a mesma regra já aplicada a `B` e `F`, só que estendida. Nenhum endpoint de PXR pode
  ser chamado de "ligação" ou "funcional" pelo sufixo.
- **Só há modo agonista na triagem.** O programa Tox21 não tem contra-ensaio de antagonismo de PXR, então
  a `direcao` do PXR nunca sai como antagonista. A máscara de citotoxicidade **não** é aplicada, porque o
  projeto só a aplica aos modos antagonistas; as colunas `citotoxico` e `citotox_confunde` são calculadas
  e gravadas assim mesmo, e 774 dos 2.759 ativos (28%) coincidem com perda de viabilidade.
- **O BindingDB não acrescenta camada de PXR.** São 2.036 registros, mas 1.530 são reimportação do ChEMBL
  e, depois do funil, sobram 215 medidas com apenas **31 moléculas** ausentes da camada do ChEMBL. A
  fração de literatura, a única que a medição do projeto mostrou melhorar o balanço de classes, são 14
  registros. Abaixo do mínimo de 50 do funil, então nenhum dataset foi gerado.

## Arquitetura: três camadas independentes

As camadas nunca são fundidas, porque medem coisas diferentes e os ensaios não são comparáveis entre si.

| Camada | Pasta | Modelos | Conteúdo |
|---|---|---|---|
| ChEMBL IC50 | `models_final/` | 18 | regressão de pIC50 e classificação, por alvo e tipo de ensaio (ligação `B` ou funcional `F`) |
| ChEMBL Ki, EC50 e Kd | `models_layers/` | 40 | mesmos alvos, cada medida em seu próprio modelo. Inclui os 8 modelos de PXR (`EC50_A`, `EC50_B`, `IC50_B`, `EC50_F`), que não têm medida de Ki nem Kd aproveitável |
| Tox21 qHTS | `models_tox21/` | 12 | classificação por ensaio (química ambiental). São 13 ensaios curados e 12 com modelo: o de TR agonista tem só 27 ativos e fica de fora. O contra-ensaio de viabilidade só se aplica aos modos antagonistas, e o de autofluorescência só aos ensaios BLA |

São **70 modelos treinados**, dos quais **52 implantáveis**.

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

Requer Python 3.11 ou 3.12.

```bash
cd ~/Downloads/pipeline_edc_qsar
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Se o `python3 -m venv` falhar com `No module named 'ensurepip'` (acontece no Ubuntu, onde o módulo vem
no pacote `python3-venv`), use o `virtualenv`, que traz o `pip` embutido e não precisa de `sudo`:

```bash
python3 -m pip install --target .tools virtualenv
PYTHONPATH=.tools python3 -m virtualenv -p python3 .venv
.venv/bin/pip install -r requirements.txt
```

Depois disso a pasta `.tools/` pode ser apagada: ela só serve para criar o ambiente.

Conferência da instalação:

```bash
# deve imprimir 2180
.venv/bin/python -c "import features as F; print(2048 + len(F.DESC_NAMES))"

# deve imprimir 'carregados 70'
.venv/bin/python -c "
import glob, joblib
print('carregados', sum(bool(joblib.load(f)) for f in glob.glob('models_*/*.joblib')))"

# previsão de referência: o alerta de ERalfa deve sair 'alto'
.venv/bin/python edc.py testar --smiles 'CC(C)(c1ccc(O)cc1)c1ccc(O)cc1'
```

As versões de `rdkit` e `scikit-learn` precisam ser as de `requirements.txt`. Com outra versão do RDKit o
número de descritores muda e o preditor falha com `Feature shape mismatch`. O `requirements.txt` explica a
situação do `xgboost`, cuja versão de treino (3.4.2) só existe no conda-forge.

### Atalho

O executável `edc` na raiz chama o `edc.py` com o ambiente da pasta, de qualquer diretório:

```bash
./edc ajuda
~/Downloads/pipeline_edc_qsar/edc testar --arquivo moleculas.csv --coluna-nome nome
```

## Uso

O ponto de entrada é o `edc.py`, com três comandos:

```bash
python edc.py ajuda                                              # lista os comandos
python edc.py explicar                                           # como usar e como interpretar
python edc.py testar --smiles "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1"    # uma molécula
python edc.py testar --smiles "SMILES_1" "SMILES_2"              # várias moléculas
python edc.py testar --arquivo moleculas.csv --coluna-nome nome  # uma lista
python edc.py testar --arquivo moleculas.csv --eixos             # com a leitura em dois eixos
```

Sem argumentos, `python edc.py` mostra a lista de comandos. Com `--saida PREFIXO` você escolhe onde
gravar; sem ele, a saída vai para `resultados/AAAAMMDD_HHMMSS_predicoes/`. O comando `explicar` também
conta quantos modelos estão instalados em cada camada.

Para a execução arquivada, que além da previsão gera os gráficos:

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

- Divisão externa por scaffold de Bemis-Murcko (80/20) e seleção de modelo por validação cruzada 5-fold no treino. **Atenção:** o código preenche o *teste* com os **maiores** grupos de scaffold (o docstring de `features.scaffold_split` afirma o contrário). O teste externo é, portanto, extrapolação para as famílias químicas dominantes, concentrada em poucos scaffolds (de 2 a 15 nos endpoints de Ki/Kd/EC50), e as métricas de teste são um limite inferior de alta variância. A validação cruzada do treino, ao contrário, é aleatória e não aninhada, então o `Q2_cv` é otimista.
- Modelo final reajustado em 100% dos dados curados, com as métricas vindas do teste externo por scaffold. As métricas gravadas dentro do artefato de produção (`Q2_cv`, `R2_teste_scaffold`) são, portanto, de **outro** ajuste: o de 80%. O limiar de domínio gravado no artefato também é recalculado sobre 100% dos dados, então é mais largo que o limiar validado.
- Y-scrambling e domínio de aplicabilidade por similaridade de Tanimoto.
- Predição conformal indutiva com cobertura nominal de 80%, e limiar de decisão otimizado por MCC em vez de 0,5. Duas ressalvas medidas: a cobertura **empírica** fica abaixo de 0,80 em 8 dos 9 endpoints de IC50; e o quantil conformal é calibrado num modelo ajustado em 75% do treino, mas aplicado ao modelo de produção ajustado em 100%, de modo que a garantia de cobertura não vale para o modelo implantado.
- Curadoria: estruturas válidas, SMILES padronizados, remoção de sais e solventes e deduplicação por InChIKey, com a mediana das réplicas e descarte do composto quando elas variam mais de 1 unidade log. Tipos de medida (IC50, Ki, Kd, EC50) e tipos de ensaio nunca são misturados. A estereoquímica é apenas **registrada** (coluna `undefined_stereo`, que conta centros sem atribuição) e nunca filtra. A deduplicação usa o InChIKey completo, que inclui a camada estereoquímica, então estereoisômeros definidos permanecem como compostos distintos, e um composto com estereoquímica indefinida convive como linha separada dos seus enantiômeros resolvidos.

## Problemas conhecidos

Verificados no código e nos dados. O detalhamento, com os números e os arquivos de diagnóstico, está em
`documentos/verificacao_auditoria_codigo.md`.

### Corrigidos

Nenhuma correção altera a previsão de modelo nenhum: o que muda é quais modelos votam, como
os votos são agregados e o que a saída informa.

**Critério de implantação, agora auditável** — `criterio_implantacao.py`

O flag `nao_implantavel` gravado nos `.joblib` de `models_final/` vinha de um script ausente do
repositório, e não era auditável: `ESR2_B` passava com R² de 0,098 e `AR_F` era reprovado com 0,076,
sem regra escrita em lugar nenhum. O critério foi reconstruído a partir das métricas que os próprios
artefatos guardam:

| Condição | Veredito |
|---|---|
| `n_train < 200` | reprovado, em qualquer camada |
| regressão com `Q2_cv < 0,30` | reprovado |
| regressão com `R2_teste_scaffold < 0` | reprovado (pior que prever a média) |
| classificação com `ROC_AUC_teste < 0,60` | reprovado |

Ele reproduz todos os flags gravados menos dois, divergindo só em `AR_F` (n = 618, Q² = 0,450,
R² = 0,076, ROC-AUC = 0,772), que passa a ser implantado. Resultado: **52 implantados de 70 treinados**,
já com os 9 modelos de PXR (5 implantáveis: `PXR_EC50_A` e `PXR_EC50_B` nas duas tarefas, mais o ensaio
Tox21; `PXR_IC50_B` com n = 176 e `PXR_EC50_F` com n = 62 são reprovados por tamanho).
`criterio_implantacao.csv` traz o motivo de cada reprovação.

Ressalva importante: **os 11 modelos do Tox21 contam como implantados por ausência de critério**, não
por terem passado em um — nem os `.joblib` nem o `LEIAME_modelos_tox21.csv` guardam métrica de
validação. Isso precisa ser dito em qualquer relatório que use a camada de triagem.

**Nova regra de níveis, baseada em unidades de evidência**

Um modelo não é uma evidência. Duas redundâncias foram medidas: regressão e classificação do mesmo
endpoint usam os mesmos compostos e as mesmas features, diferindo só pela binarização em pAct ≥ 6; e
os ensaios Tox21 de um mesmo receptor e modo são a mesma biblioteca (até 99,8% em comum) em linhagens
celulares diferentes. Os votos passam a ser colapsados antes de qualquer contagem:

- **eixo farmacológico:** 1 unidade por endpoint do ChEMBL (receptor × medida × tipo de ensaio). A
  classificação decide a positividade; a regressão entra como potência, não como voto separado.
- **eixo de triagem:** 1 unidade por receptor e modo, colapsando linhagens.

A escada passou a ser:

| Nível | Condição |
|---|---|
| `nao avaliado` | nenhuma unidade avaliada, ou uma só e negativa |
| `sem alerta` | todas as unidades avaliadas negativas (2 ou mais) |
| `alto` | os **dois** eixos positivos (convergência) |
| `moderado` | ≥ 2 unidades independentes positivas num eixo só |
| `baixo` | exatamente 1 unidade positiva |

O estado `nao avaliado` é o ganho principal: a regra anterior chamava de "sem alerta" tanto a molécula
testada e negativa quanto aquela que modelo nenhum conseguiu avaliar. No `painel_ref.csv`, **14 dos 40
pares que eram "sem alerta" eram, na verdade, não avaliados**. Os 10 "alto" permaneceram "alto", todos
por convergência.

**Saída**

- O `_detalhado.csv` traz **apenas os modelos que votaram**. `--detalhe implantados` acrescenta os que
  não votaram por estarem fora do domínio, e `--detalhe todos` acrescenta também os reprovados. Para o
  bisfenol A isso dá 27, 47 e 61 linhas. A tabela completa continua sendo usada internamente pelo
  `--eixos`, cuja cobertura depende de comparar as unidades implantáveis com as avaliadas.
- Coluna nova `motivo_nao_voto`: `modelo reprovado (<motivo>)` ou `molecula fora do dominio`.
- Toda combinação molécula × receptor ganha linha no `_resumo.csv`, inclusive as sem nenhuma cobertura,
  que saem como `nao avaliado`. Antes, um receptor cujas unidades ficassem todas fora do domínio
  desaparecia da saída.
- A coluna `evidencia` fala só do que foi avaliado (`positivas 1/3`). A cobertura continua no CSV, nas
  colunas `unidades_implantaveis` e `cobertura`, mas sai da tabela impressa — a decisão foi deixar a
  tela com o que de fato foi medido. Vale consultá-las antes de comparar duas moléculas: no painel, 33
  dos 83 pares avaliados pertencem a grupos de mesmo nível e mesma fração positiva cuja cobertura varia
  mais de 0,40 (`Ibuprofeno · ERβ` tem cobertura 0,25 e `17β-Estradiol · TR` tem 1,00, e ambos saem
  "positivas 1/1"). A linha impressa acima da tabela dá a cobertura agregada da execução.
- O `_resumo.csv` passou a `evidencia`, `unidades_positivas`, `unidades_avaliadas`,
  `unidades_implantaveis`, `fracao_positiva`, `cobertura`, `concordancia`, `direcao` e `alerta`.
- A tela informa quantas unidades de evidência puderam votar, do total disponível.

**Produziam resultado errado sem avisar**

- **Camada identificada pelo nome da pasta.** A comparação passou a usar o nome *base*
  (`os.path.basename`), então `--model-dirs` com caminho absoluto funciona, e uma pasta com nome
  desconhecido levanta erro. Antes, `--model-dirs /abs/.../models_tox21` rotulava os 11 modelos do
  Tox21 como `ChEMBL AR/ESR1/PPARG/TR`, com o modo do ensaio virando nome de receptor.
- **Nomes repetidos na entrada.** O `predict_edc2.py` agora levanta erro explícito, em vez de fundir
  moléculas diferentes numa linha do resumo.
- **Moléculas sem scaffold de Bemis-Murcko.** Cada uma passou a receber o próprio grupo. O teste
  externo de `TOX21_AR_antagonista_HEK293T` caiu de 95,2% para 0,0% de acíclicas e o de
  `TOX21_TR_antagonista_GH3` de 91,9% para 0,1%; `PPARG_Kd_B`, cujo teste tinha uma **única classe**,
  passou a ter duas. Afeta apenas retreino.
- **Queda silenciosa para divisão aleatória.** `features.scaffold_split` levanta erro, e só aceita
  divisão aleatória com `permitir_aleatorio=True`.

**Impediam rodar de ponta a ponta**

- **BindingDB.** O `pipeline_curate.py` lê o `.parquet` que o `extract_bindingdb.py` grava e
  **constrói** o `cache/bindingdb_std_lookup.pkl`, em vez de ler um arquivo que nenhum script criava.
- **Nomes das pastas de dados.** `modeling_layers.py` e `modeling_tox21.py` aceitam tanto `data_layers/`
  e `data_tox21/` quanto `datasets_camadas_por_endpoint/` e `datasets_tox21_por_ensaio/`.

### Em aberto

**Reprodutibilidade do que está implantado** — não é corrigível por mudança de código:

- Os 18 modelos de `models_final/` **não são regeráveis por nenhum script da pasta**. O
  `LEIAME_calibracao.csv` tem colunas (`limiar_producao`, `ganho_limiar`, `conformal_modo`,
  `cobertura_scaffold`) que o `conformal.py` não produz; ele grava outro arquivo,
  `results/conformal_e_limiares.csv`. O script que fecha essa etapa não está no repositório, e reescrevê-lo
  produziria modelos equivalentes, não os mesmos. (A quinta coluna, `nao_implantavel`, deixou de depender
  desse script: o `criterio_implantacao.py` a reconstrói com regra escrita.)
- `finalize_models.py` e `conformal.py` dependem de `data/dataset_*.csv` (datasets curados de IC50) e de
  `results/desempenho_modelos.csv`, ambos ausentes. Os datasets de IC50 ficaram deliberadamente fora do
  repositório.

**Limites estruturais do consenso** — dependem de retreino ou de dados novos, não de código:

- Os votos continuam sem peso: duas unidades positivas valem o mesmo, venham de 2.700 compostos
  (`ESR1_B`) ou de 232 (`AR_EC50_B`). Ponderar exigiria uma medida de confiança comparável entre as
  camadas, que a de triagem não tem.
- O teto do alerta depende do receptor: "alto" só é atingível em AR, ERα, PPARγ e PXR, porque exige um
  voto de cada camada. ERβ não tem ensaio Tox21 e chega no máximo a "moderado"; THRβ chega a "baixo";
  THRα não tem nenhum modelo ChEMBL implantável e não aparece no resumo. O PXR é o mais magro dos que
  alcançam "alto": 2 unidades farmacológicas e 1 de triagem, contra 4 e 4 do ERα. Um "sem alerta" não significa o mesmo em
  receptores diferentes. A coluna `modelos_implantaveis_receptor` deixa esse teto visível.
- **No PXR, os dois eixos não estão em pé de igualdade no painel de referência.** Das 18 moléculas do
  `painel_ref.csv`, **14 estão no conjunto de treino do Tox21 PXR e nenhuma no do ChEMBL PXR** (Tanimoto
  mediano ao vizinho de treino: 1,00 na triagem contra 0,28 no eixo farmacológico). Por isso 13 das 18
  aparecem como "divergente: triagem positiva, farmacológico negativo" — o modelo de triagem reconhece
  compostos que já viu, e os farmacológicos estão extrapolando. Esse padrão é da composição do painel, não
  uma propriedade do PXR, e o painel não serve para estimar desempenho do PXR em moléculas novas.
- "Ativo" não quer dizer o mesmo em cada camada: no ChEMBL é potência de até 1 µM com limiar de
  probabilidade 0,5; no Tox21 é a chamada do depositante, com limiares otimizados de 0,11 a 0,40. Por
  isso os dois eixos nunca são somados numa fração única.
- A camada Tox21 não tem métrica de validação gravada, então seus 11 modelos nunca foram submetidos ao
  critério de implantação. São 11 das 47 unidades que votam.
- O `painel_ref.csv` **não é um conjunto de avaliação cego**: 59,9% dos votos que contam vêm de modelos
  que já tinham a molécula no próprio treino (84,3% no eixo de triagem). Ele serve para conferir o
  mecanismo do consenso, nunca para estimar acerto em moléculas novas.

O `consenso_v2/` reescreve só a agregação e deixa esses pontos explícitos no relatório, sem retreinar
nada. Veja `documentos/consenso_v2_especificacao.md`.

## Limitações

- Triagem computacional. Um alerta indica plausibilidade e não substitui ensaio experimental nem avaliação regulatória.
- Moléculas fora do domínio de aplicabilidade não entram na regra de alerta, mesmo quando o modelo devolve uma probabilidade alta. Compostos de classes pouco representadas no treino, como ftalatos, podem ficar sem cobertura.
- Mecanismos que não passam pelos receptores modelados não são detectados. Por exemplo, efeitos via aromatase (CYP19A1) ou esteroidogênese exigiriam endpoints adicionais.
- A camada Tox21 é de química ambiental e usa outro tipo de medida (AC50) que as camadas do ChEMBL.

## Fontes de dados

ChEMBL 37, BindingDB e Tox21. Os dados derivados seguem os termos de cada fonte. O ChEMBL é distribuído sob CC BY-SA 3.0. Consulte os termos das demais fontes antes de redistribuir os datasets.

## Licença

A definir.
