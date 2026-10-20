# Modelos QSAR para desregulação endócrina — AR, ERα, ERβ, THRα, THRβ e PPARγ

Dataset curado de IC50 e modelos de aprendizado de máquina para triagem de novos
compostos com potencial de interação com receptores nucleares associados à
desregulação endócrina.

---

## 1. Fontes de dados

| Fonte | Situação | Registros de IC50 nos 6 alvos |
|---|---|---|
| **ChEMBL 37** (dump SQLite oficial, SHA256 conferido) | fonte primária | 19.472 |
| **BindingDB 2026-09** (dump TSV completo) | fonte secundária / validação externa | 17.544 |
| **PubChem BioAssay** (SDQ, por accession UniProt) | consultado, sem contribuição independente | 20.882, 100% depositados pelo ChEMBL ou BindingDB |
| **Cortellis** (Clarivate) | **não acessível** | — |

**Sobre o Cortellis.** É um produto comercial licenciado, sem API pública e sem
credencial disponível neste ambiente. Ficou fora da coleta. Caso haja acesso
institucional, uma exportação em CSV pode ser integrada à mesma estrutura de
proveniência sem alteração do pipeline.

**Sobre a API REST do ChEMBL.** Durante a execução, os endpoints de lista
responderam HTTP 500 em cerca de 70% das chamadas. A extração foi migrada para
o dump SQLite oficial (`chembl_37_sqlite.tar.gz`, 5,76 GB), consultado por SQL —
solução mais completa e reprodutível do que a paginação da API.

**Alvos.** Resolvidos por accession UniProt, restritos a `SINGLE PROTEIN` humano:

| Receptor | ChEMBL ID | UniProt |
|---|---|---|
| Receptor de andrógeno (AR) | CHEMBL1871 | P10275 |
| Receptor de estrógeno α (ERα/ESR1) | CHEMBL206 | P03372 |
| Receptor de estrógeno β (ERβ/ESR2) | CHEMBL242 | Q92731 |
| PPAR-γ | CHEMBL235 | P37231 |
| Receptor de hormônio tireoidiano α (THRα) | CHEMBL1860 | P10827 |
| Receptor de hormônio tireoidiano β (THRβ) | CHEMBL1947 | P10828 |

Complexos proteicos, famílias e grupos de seletividade (ex.: "Estrogen receptor"
como família, "PPAR alpha/gamma") foram **excluídos** — o alvo molecular precisa
ser inequívoco para que o endpoint seja interpretável.

---

## 2. Comparabilidade dos ensaios

O requisito central foi não misturar medidas biologicamente incomparáveis.
Filtros aplicados, em ordem, sobre os registros do ChEMBL:

| filtro                                                  |   registros restantes |
|:--------------------------------------------------------|----------------------:|
| 0. IC50 brutos (alvo humano SINGLE PROTEIN)             |                 19472 |
| 1. alvo Homo sapiens                                    |                 19472 |
| 2. confidence_score >= 8 (alvo proteico direto)         |                 19472 |
| 3. standard_relation '='                                |                 14575 |
| 4. unidade conversível para nM                          |                 14572 |
| 5. valor em faixa plausível (pIC50 2-12)                |                 14568 |
| 6. sem data_validity_comment problemático               |                 14513 |
| 7. com estrutura química definida                       |                 14502 |
| 8. estrutura padronizada válida (sal/solvente removido) |                 14445 |

Além disso:

- **Separação por tipo de ensaio.** Ensaios de **ligação** (`assay_type = B`,
  tipicamente deslocamento de radioligante) e **funcionais** (`assay_type = F`,
  tipicamente gene repórter em célula) geram datasets **independentes**. Um IC50
  de deslocamento competitivo e um IC50 de inibição transcricional medem
  fenômenos distintos e não foram fundidos em nenhum momento.
- **confidence_score ≥ 8** — atribuição direta a proteína única. Todos os
  19.472 registros já atendiam a esse critério.
- **Apenas relação `=`.** Registros censurados (`>`, `<`) foram descartados:
  4.897 registros, 25% do total bruto.
- **Modo de ação anotado** (agonista/antagonista/indefinido) a partir da
  descrição do ensaio, disponível como coluna para filtragem adicional. Os
  conjuntos funcionais são dominados por um dos modos (AR funcional: 98,7%
  antagonista; PPARγ funcional: 7,5% antagonista, ou seja, predominantemente
  agonista) — misturar esses dois receptores num mesmo modelo funcional seria
  incorreto, e eles permanecem separados.
- **Réplicas.** Mediana por InChIKey dentro de cada endpoint; compostos com
  amplitude entre réplicas superior a **1 unidade log** foram descartados
  (747 compostos no total, sendo 518 em ERα ligação — 16% daquele conjunto).

### Padronização estrutural (RDKit)

`Cleanup` → seleção do maior fragmento (remoção de sais e solventes) →
normalização de grupos funcionais → neutralização de cargas → canonicalização
tautomérica → InChIKey como identificador de identidade. Rejeitados: estruturas
inválidas, inorgânicos, elementos fora do conjunto orgânico, menos de 6 átomos
pesados e massa molecular fora de 100–1000 Da (30 estruturas no ChEMBL, 66 no
BindingDB). Estereocentros indefinidos foram **anotados** (coluna
`undefined_stereo`), não removidos: 2.153 dos 9.110 pares composto–alvo.

---

## 3. Datasets resultantes

**7245 compostos únicos, 9.110 pares composto–endpoint, 9 endpoints.**

| receptor   | ensaio    |   registros IC50 |   compostos |   descartados (réplicas discordantes) |   compostos finais |   pIC50 mediana |   fração modo antagonista |
|:-----------|:----------|-----------------:|------------:|--------------------------------------:|-------------------:|----------------:|--------------------------:|
| ERα        | ligação   |             5637 |        3218 |                                   518 |               2700 |            7    |                     0.413 |
| AR         | ligação   |             2687 |        1883 |                                    76 |               1807 |            6.46 |                     0.603 |
| ERβ        | ligação   |             2146 |        1608 |                                    60 |               1548 |            6.89 |                     0.501 |
| PPARγ      | ligação   |             1665 |        1464 |                                    27 |               1437 |            6.48 |                     0.153 |
| AR         | funcional |              888 |         656 |                                    38 |                618 |            6.7  |                     0.987 |
| THRβ       | ligação   |              519 |         402 |                                     6 |                396 |            6.14 |                     0.831 |
| PPARγ      | funcional |              300 |         231 |                                     3 |                228 |            6.28 |                     0.075 |
| THRα       | ligação   |              272 |         221 |                                     1 |                220 |            5.8  |                     0.714 |
| ERα        | funcional |              220 |         174 |                                    18 |                156 |            6.38 |                     0.788 |

![Visão geral dos datasets curados]({{artifact:art_f6165da8-e2cb-44de-879e-aa8e7ecf2813}})

*(a) Tamanho dos conjuntos. (b) Distribuição de pIC50 com o limiar de atividade
de 1 µM. (c) t-SNE de ECFP4 dos conjuntos de ligação — cada receptor ocupa uma
região química distinta. (d) Compostos testados em pares de receptores: a
sobreposição se concentra em ERα/ERβ (1.170 compostos) e THRα/THRβ (200), como
esperado pela homologia dos pares.*

Cada linha do CSV carrega proveniência completa: `molecule_chembl_id`, nome,
SMILES original e padronizado, InChIKey, tipo/valor/unidade da medida, pIC50,
pChEMBL, `assay_chembl_id` e descrição textual do ensaio, `confidence_score`,
BAO format, organismo, `target_chembl_id`, documento de origem, DOI, PMID, ano e
periódico. 83% dos compostos têm DOI e PMID rastreáveis; as colunas `all_assays`,
`all_docs`, `all_dois` e `all_pmids` listam **todas** as fontes que sustentam o
valor agregado de cada composto.

---

## 4. Modelagem

**Representação.** ECFP4 contado (raio 2, 2048 bits) concatenado a 201
descritores 2D do RDKit. Para modelos lineares e de margem, o pipeline inclui
`VarianceThreshold` e padronização — ajustados apenas no treino, dentro da
validação cruzada.

**Protocolo.**

1. Split externo por **scaffold de Bemis-Murcko** (80/20): o conjunto de teste
   contém apenas esqueletos moleculares ausentes do treino. É um teste
   deliberadamente severo, próximo do cenário real de priorizar química nova.
2. Seleção de modelo por busca em grade com **CV 5-fold** dentro do treino
   (Ridge, SVR-RBF, Random Forest, XGBoost para regressão; Regressão Logística,
   SVC-RBF, Random Forest, XGBoost para classificação).
3. Refit no treino completo e avaliação no teste externo.
4. **y-scrambling** (5–10 permutações) sobre o melhor modelo.
5. **Domínio de aplicabilidade**: maior similaridade de Tanimoto (ECFP4) ao
   conjunto de treino; o limiar é o percentil 5 da distribuição de similaridade
   ao vizinho mais próximo dentro do treino.
6. **Modelos de produção** reajustados em 100% dos dados curados. As métricas
   reportadas continuam vindo do protocolo com teste externo.

**Classificação.** Limiar de atividade em pIC50 ≥ 6 (IC50 ≤ 1 µM), convenção
usual em triagem de desreguladores endócrinos. A fração de ativos varia de 0,41
(THRα) a 0,79 (AR funcional).

| receptor   | ensaio    |   n_treino_final | modelo_regressao   |   Q2_cv |   R2_teste_scaffold | modelo_classificacao   |   AD_limiar_tanimoto |
|:-----------|:----------|-----------------:|:-------------------|--------:|--------------------:|:-----------------------|---------------------:|
| AR         | ligação   |             1807 | XGBoost            |   0.604 |               0.273 | XGBoost                |                0.5   |
| AR         | funcional |              618 | RandomForest       |   0.45  |               0.076 | RandomForest           |                0.537 |
| ERα        | ligação   |             2700 | XGBoost            |   0.793 |               0.66  | XGBoost                |                0.327 |
| ERα        | funcional |              156 | RandomForest       |   0.589 |              -0.257 | RandomForest           |                0.324 |
| ERβ        | ligação   |             1548 | XGBoost            |   0.602 |               0.098 | RandomForest           |                0.297 |
| PPARγ      | ligação   |             1437 | XGBoost            |   0.662 |               0.477 | Logistic               |                0.555 |
| PPARγ      | funcional |              228 | Ridge              |   0.627 |               0.376 | SVC                    |                0.494 |
| THRα       | ligação   |              220 | RandomForest       |   0.72  |              -0.065 | RandomForest           |                0.571 |
| THRβ       | ligação   |              396 | XGBoost            |   0.767 |               0.58  | XGBoost                |                0.619 |

---

## 5. Desempenho

![Validação dos modelos]({{artifact:art_69bd504f-ddad-4ae1-985f-cce7a0ed11e0}})

### Regressão de pIC50

| receptor   | ensaio    | modelo       |   Q2_cv |   RMSE_cv |   R2_teste |   RMSE_teste |   MAE_teste |
|:-----------|:----------|:-------------|--------:|----------:|-----------:|-------------:|------------:|
| AR         | ligação   | XGBoost      |   0.604 |     0.632 |      0.273 |        0.776 |       0.604 |
| AR         | funcional | RandomForest |   0.45  |     0.634 |      0.076 |        0.662 |       0.523 |
| ERα        | ligação   | XGBoost      |   0.793 |     0.717 |      0.66  |        0.747 |       0.554 |
| ERα        | funcional | RandomForest |   0.589 |     0.759 |     -0.257 |        1.068 |       0.867 |
| ERβ        | ligação   | XGBoost      |   0.602 |     0.771 |      0.098 |        0.961 |       0.779 |
| PPARγ      | ligação   | XGBoost      |   0.662 |     0.671 |      0.477 |        0.777 |       0.623 |
| PPARγ      | funcional | Ridge        |   0.627 |     0.625 |      0.376 |        0.947 |       0.756 |
| THRα       | ligação   | RandomForest |   0.72  |     0.535 |     -0.065 |        1.298 |       0.91  |
| THRβ       | ligação   | XGBoost      |   0.767 |     0.697 |      0.58  |        0.427 |       0.339 |

### Classificação ativo/inativo (pIC50 ≥ 6)

| receptor   | ensaio    | modelo       |   ROC_AUC_cv |   ROC_AUC_teste |   PR_AUC_teste |   BalAcc_teste |   MCC_teste |
|:-----------|:----------|:-------------|-------------:|----------------:|---------------:|---------------:|------------:|
| AR         | ligação   | XGBoost      |        0.88  |           0.759 |          0.851 |          0.693 |       0.391 |
| AR         | funcional | RandomForest |        0.844 |           0.772 |          0.984 |          0.5   |       0     |
| ERα        | ligação   | XGBoost      |        0.951 |           0.855 |          0.965 |          0.749 |       0.485 |
| ERα        | funcional | RandomForest |        0.872 |           0.94  |          0.974 |          0.875 |       0.831 |
| ERβ        | ligação   | RandomForest |        0.894 |           0.821 |          0.965 |          0.748 |       0.391 |
| PPARγ      | ligação   | Logistic     |        0.904 |           0.899 |          0.949 |          0.815 |       0.706 |
| PPARγ      | funcional | SVC          |        0.909 |           0.932 |          0.949 |          0.91  |       0.821 |
| THRα       | ligação   | RandomForest |        0.957 |           0.528 |          0.811 |          0.585 |       0.284 |
| THRβ       | ligação   | XGBoost      |        0.918 |           0.841 |          0.736 |          0.783 |       0.486 |

### Validação externa independente — BindingDB

817 compostos presentes no BindingDB e **ausentes** do ChEMBL, ou seja, nunca
vistos pelos modelos em nenhuma etapa:

| receptor   |   n_externo |   frac_no_dominio |     R2 |   RMSE |   ROC_AUC |   BalAcc |   MCC |
|:-----------|------------:|------------------:|-------:|-------:|----------:|---------:|------:|
| ERα        |         378 |             0.944 |  0.502 |  1.272 |     0.883 |    0.703 | 0.56  |
| ERβ        |         235 |             0.702 | -0.286 |  1.286 |     0.761 |    0.77  | 0.519 |
| AR         |         128 |             0.68  |  0.466 |  0.917 |     0.837 |    0.811 | 0.682 |
| PPARγ      |          44 |             0.682 |  0.148 |  0.799 |     0.841 |    0.736 | 0.388 |
| THRα       |          17 |             0.824 |  0.888 |  0.359 |     1     |    1     | 1     |
| THRβ       |          15 |             0.133 |  0.362 |  0.689 |     1     |    0.964 | 0.681 |

### Leitura dos resultados

- **O y-scrambling exclui correlação por acaso** em todos os endpoints: Q²
  entre −0,43 e −0,14 com rótulos permutados, contra 0,45 a 0,79 com os rótulos
  reais.
- **A classificação generaliza melhor que a regressão.** Fora do treino, o
  ROC-AUC fica em 0,76–0,94 em 8 dos 9 endpoints, enquanto o R² de regressão cai
  de 0,45–0,79 (CV) para −0,26 a 0,66 (scaffold). Para priorização de compostos
  novos, use a **probabilidade de atividade**; trate o pIC50 previsto como
  ordenação aproximada, não como estimativa quantitativa.
- **Quatro modelos são sólidos**: ERα ligação (R² 0,66 no teste por scaffold,
  RMSE 0,75 log), THRβ ligação (0,58 / 0,43), PPARγ ligação (0,48 / 0,78) e AR
  ligação (0,27 / 0,78 — fraco em regressão, ROC-AUC 0,76 em classificação).
- **Três modelos são frágeis e devem ser usados só como alerta qualitativo**:
  ERα funcional (n = 156), THRα ligação (n = 220, ROC-AUC 0,53 no teste por
  scaffold) e AR funcional (acurácia balanceada 0,50 — o classificador prevê a
  classe majoritária, que representa 79% do conjunto).
- **THRα e THRβ têm dados demais em comum e de menos no total.** 200 dos 220
  compostos de THRα também estão em THRβ; os modelos dos dois receptores não são
  independentes o suficiente para sustentar previsão de seletividade α/β.
- **O domínio de aplicabilidade separa o erro de forma clara apenas em THRα**
  (RMSE 1,01 dentro contra 1,98 fora). Nos demais endpoints o limiar do
  percentil 5 é permissivo demais para funcionar como filtro de confiança —
  interprete a similaridade de Tanimoto reportada como valor contínuo, não como
  binário.

---

## 6. Uso do preditor

```bash
python predict_edc.py --smiles "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1" --model-dir models_final
python predict_edc.py --input moleculas.csv --smiles-col smiles --name-col nome \
                      --model-dir models_final --out predicoes.csv
```

Para cada molécula e cada endpoint a saída traz `pIC50_prev`, `prob_ativo`,
`tanimoto_viz` (similaridade ao vizinho mais próximo do treino), `no_dominio`, e
um resumo multi-alvo (`n_endpoints_ativos_no_dominio`, `perfil_multialvo`).

**Demonstração com desreguladores endócrinos conhecidos** (arquivo
`predicoes_demo.csv`) ilustra tanto o alcance quanto o limite do modelo. O
bisfenol A é recuperado como ativo em ERα e ERβ (pIC50 previsto 5,96 e 6,02,
Tanimoto 1,0 — ele está nos conjuntos de treino) e a flutamida como ativa em AR.
Já a genisteína é subestimada em ERα (4,94 previsto contra ≈7,5 na literatura) e
a rosiglitazona em PPARγ (4,75 contra ≈8), **porque nenhuma das duas está nos
conjuntos**: a potência delas é reportada como EC50, Ki ou afinidade relativa, e
o recorte IC50 as exclui. DDT, atrazina e DEHP caem fora do domínio de
aplicabilidade em quase todos os endpoints.

---

## 7. Limitações

1. **O recorte IC50 é farmacológico, não toxicológico.** A química do ChEMBL é
   dominada por séries de descoberta de fármacos. Contaminantes ambientais
   clássicos (ftalatos, organoclorados, triazinas) estão mal representados e
   caem fora do domínio de aplicabilidade. Para triagem ambiental, o caminho
   natural é complementar com os ensaios qHTS do Tox21/ToxCast — que reportam
   **AC50**, um endpoint distinto, e por isso exigiriam um modelo separado. Essa expansão é detalhada na seção 9.
2. **Compostos ausentes por endpoint.** Ligantes bem conhecidos ficam de fora
   quando sua potência foi publicada como EC50/Ki em vez de IC50 (genisteína,
   rosiglitazona). Um dataset paralelo de Ki e outro de EC50 ampliariam a
   cobertura — sem nunca fundi-los com o de IC50.
3. **Variabilidade experimental herdada.** Mesmo após o corte de réplicas
   discordantes, a dispersão entre laboratórios impõe um piso de erro da ordem
   de algumas décimas de unidade log.
4. **Agonismo e antagonismo não são distinguidos na predição.** O modo de ação
   está anotado no dataset, mas os modelos preveem potência, não direção do
   efeito — um ponto relevante, já que a consequência toxicológica de agonizar ou
   antagonizar um receptor nuclear é oposta.
5. **Ausência do Cortellis** e de fontes de literatura não indexadas em ChEMBL
   ou BindingDB.
6. **Uso pretendido.** Priorização de compostos para ensaio experimental. Não
   substitui teste in vitro nem constitui avaliação regulatória de risco.

---

## 8. Reprodutibilidade

| Arquivo | Conteúdo |
|---|---|
| `extract_chembl.py` | SQL de extração do dump SQLite do ChEMBL 37 |
| `extract_bindingdb.py` | filtro do dump completo do BindingDB por accession |
| `curation.py` | padronização RDKit, conversão para pIC50, agregação de réplicas |
| `pipeline_curate.py` | pipeline de curadoria e montagem dos datasets |
| `features.py` | ECFP4 + descritores, split por scaffold, Tanimoto |
| `modeling.py` | treino, seleção, validação, y-scrambling, domínio de aplicabilidade |
| `finalize_models.py` | refit dos modelos de produção em 100% dos dados |
| `predict_edc.py` | predição multi-alvo para novos SMILES |

Ambiente: Python 3.11, RDKit 2026.03.6, scikit-learn 1.9.1, XGBoost, pandas,
numpy. Semente fixa (42) em todos os splits e modelos.

---

## 9. Expansão de escopo: modelo independente para AC50 (Tox21/ToxCast)

Esta seção responde à limitação 7.1 e define o escopo de uma nova frente de modelagem. Ela é uma
proposta: nenhum resultado abaixo foi obtido ainda.

### 9.1 Motivação

O recorte IC50 do ChEMBL é dominado por séries de descoberta de fármacos, e contaminantes ambientais
(ftalatos, organoclorados, triazinas) ficam fora do domínio de aplicabilidade. Os ensaios de triagem de
larga escala em concentração-resposta (qHTS) do Tox21 e do ToxCast cobrem justamente a química
ambiental e industrial. Eles medem a **AC50**, a concentração em que a resposta atinge 50% da sua
amplitude máxima.

### 9.2 Estado atual no projeto

A camada Tox21 existente (`pipeline_tox21.py`, `modeling_tox21.py`, `models_tox21/`) contém
classificadores ativo/inativo por ensaio. A AC50 é convertida em pAC50 (6 − log10 AC50 em µM) e só é
usada quando o desfecho é "Active". Ainda não existe modelo que **preveja o valor** da pAC50. O
ToxCast (invitrodb) também não foi incorporado. A expansão proposta cobre essas duas lacunas.

### 9.3 Por que a AC50 não pode ser fundida ao modelo de IC50

1. **Natureza do parâmetro.** A IC50 do ChEMBL descreve inibição em ensaios curados de literatura, de
   ligação ou funcionais, com protocolos heterogêneos. A AC50 é o parâmetro de um ajuste de curva
   (Hill) sobre uma grade fixa de concentrações, gerado por um pipeline uniforme de triagem. Os
   dois números não estão na mesma escala de significado, mesmo quando as unidades coincidem.
2. **Direção do efeito.** A AC50 do qHTS é reportada para ativação (modo agonista) e para inibição
   (modo antagonista). Uma mesma AC50 pode significar efeitos opostos, e o modelo de IC50 não
   distingue direção (limitação 7.4).
3. **Censura dos dados.** Substâncias inativas não têm AC50 bem definida, ou têm um valor limite
   artificial. A distribuição de pAC50 observada é truncada pela faixa de concentração testada.
4. **Confundimento por citotoxicidade.** No modo antagonista, a queda de sinal pode refletir morte
   celular. O Tox21 traz contra-ensaios de viabilidade, que a camada IC50 não tem.
5. **Tecnologia e linhagem.** Repórter de luciferase ou beta-lactamase, linhagem celular e presença de
   autofluorescência introduzem viés específico de ensaio, ausente dos dados de ligação.
6. **Efeito do pooling.** Misturar os dois alvos de regressão deslocaria a distribuição do alvo,
   contaminaria o domínio de aplicabilidade e invalidaria a calibração conformal, que é feita por
   camada. O projeto já segue a regra de nunca fundir IC50, Ki, Kd e EC50, e a AC50 segue a mesma regra.

### 9.4 Arquitetura proposta (independente da camada IC50)

| Componente | Proposta |
|---|---|
| Unidade de modelagem | um modelo por ensaio e modo (agonista ou antagonista), sem pooling entre ensaios |
| Estrutura em dois estágios | estágio 1: classificação ativo/inativo (já existe, será reaproveitada como ponto de partida). Estágio 2: regressão de pAC50 apenas sobre os ativos |
| Censura | tratar valores limite como dado censurado (por exemplo, regressão com censura ou exclusão explícita com relato), nunca como valor observado |
| Filtros de qualidade | excluir ativos confundidos por citotoxicidade ou autofluorescência, usando os contra-ensaios já presentes em `pipeline_tox21.py` |
| Features | as mesmas de `features.py` (ECFP4 por contagem e descritores RDKit), para manter compatibilidade operacional |
| Validação | divisão externa por scaffold de Bemis-Murcko (80/20), CV 5-fold na seleção, y-scrambling, como nas outras camadas |
| Incerteza e domínio | predição conformal e domínio de aplicabilidade calibrados **somente** nos dados de AC50 |
| Saída | nova camada `models_ac50/`, com coluna própria no `_detalhado.csv`; nunca entra na média da fração ativa das camadas de IC50 sem regra explícita |

### 9.5 Dados

- **Tox21:** já usado via PubChem, com os ensaios listados em `datasets_tox21_por_ensaio/`.
- **ToxCast (invitrodb):** a incorporar. A versão exata e o mapeamento de componentes de ensaio para
  os receptores do escopo (AR, ESR1, ESR2, PPARG, THRA, THRB) devem ser confirmados na fonte antes do
  download. A comparabilidade entre ensaios será verificada no mesmo critério das demais camadas:
  mesmo alvo, mesma direção de efeito, mesma tecnologia, e nada fundido por conveniência.
- **Curadoria:** a mesma padronização de `curation.py` (sais e solventes, neutralização, tautômero,
  InChIKey), com deduplicação por InChIKey e registro de proveniência por medida.

### 9.6 Critérios de aceite

1. Desempenho no teste externo por scaffold reportado por ensaio, com intervalo e comparação contra
   um modelo de referência trivial (média do treino).
2. Cobertura de domínio medida em um painel ambiental (ftalatos, organoclorados, triazinas, bisfenóis),
   comparada com a cobertura atual das camadas de IC50. Hipótese a testar: o ganho de cobertura
   ocorre nessas classes.
3. Controles positivos: os controles positivos do painel (rosiglitazona em PPARG, genisteína em
   ESR1) devem continuar com alerta alto após a integração.
4. Teste de não fusão: verificar que a camada nova não altera o resumo de alertas existente quando
   desligada.

### 9.7 Riscos

- Poucos ativos por ensaio podem inviabilizar a regressão de alguns pares receptor-modo.
- Ensaios Tox21 e ToxCast para o mesmo alvo podem divergir, o que exige analisar a concordância
  antes de qualquer união entre as duas fontes.
- A AC50 de ativação celular não equivale a afinidade de ligação, e isso deve estar explícito na
  interpretação do alerta.
