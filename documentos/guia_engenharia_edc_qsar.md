# Engenharia do pipeline EDC-QSAR: dados, modelos e consenso

Este documento descreve o que o código em `~/Downloads/pipeline_edc_qsar` faz, na ordem em que os dados fluem. As afirmações vêm da leitura dos scripts (`curation.py`, `features.py`, `pipeline_curate.py`, `pipeline_layers.py`, `pipeline_tox21.py`, `modeling.py`, `modeling_layers.py`, `modeling_tox21.py`, `finalize_models.py`, `conformal.py`, `predict_edc2.py`), dos `LEIAME_*.csv` que acompanham os modelos, dos 61 arquivos `.joblib` carregados e de duas execuções em `resultados/`. Os datasets curados da camada IC50 não estão na pasta; para IC50 as tabelas usam os `LEIAME` e os metadados gravados nos modelos.

## 1. Visão geral

```
ChEMBL 37 (SQLite) --+-- IC50 ---------> pipeline_curate.py -> dataset_<alvo>_<B|F>.csv --> modeling.py -> finalize_models.py -> conformal.py --> models_final/  (18)
                     +-- Ki, Kd, EC50 -> pipeline_layers.py -> dataset_<alvo>_<medida>_<B|F>.csv -> modeling_layers.py --------------------------> models_layers/ (32)
Tox21 (PubChem qHTS) -> pipeline_tox21.py -> dataset_TOX21_<receptor>_<modo>_<linhagem>.csv -> modeling_tox21.py ----------------------------------> models_tox21/ (11)
BindingDB ----------> pipeline_curate.py -> externo_bindingdb.csv (validação externa dos modelos de IC50)

61 modelos --> predict_edc2.py --> 1 voto binário por modelo --> alerta por molécula x receptor
```

Regra de arquitetura: cada tipo de medida (IC50, Ki, Kd, EC50, AC50) e cada tipo de ensaio (ligação `B`, funcional `F`) tem o próprio dataset e o próprio par de modelos. Nada é fundido antes do consenso.

## 2. Dados

### 2.1 Por que as medidas não se misturam

| Medida | O que mede | Dependência do ensaio |
|---|---|---|
| IC50 | concentração que inibe 50% de uma resposta | alta: varia com a concentração de ligante/substrato e com o formato do ensaio (relaciona-se ao Ki por Cheng-Prusoff em ensaios de ligação competitiva) |
| Ki / Kd | constante de dissociação no equilíbrio | baixa em princípio (grandeza termodinâmica) |
| EC50 | concentração de efeito semimáximo em ensaio funcional | alta: depende do sistema celular, do repórter e da amplificação do sinal |
| AC50 (Tox21) | concentração de atividade semimáxima, vinda de ajuste de Hill numa grade fixa de concentrações | cobre ativação e inibição; é reportada também para inativos; no modo antagonista sofre confundimento por citotoxicidade |

As escalas numéricas coincidem: pIC50 = 9 − log10(nM) e pAC50 = 6 − log10(µM) são ambos −log10(M). Por isso um limiar "6" parece comparável entre camadas, mas o significado biológico não é o mesmo. Esse é o motivo de manter as camadas separadas e de só integrá-las no consenso.

### 2.2 Curadoria ChEMBL (IC50 em `pipeline_curate.py`; Ki, Kd e EC50 em `pipeline_layers.py`)

Os dois scripts aplicam o mesmo funil. A única diferença é o `standard_type` consultado.

1. Alvos: AR (P10275), ERα (P03372), ERβ (Q92731), PPARγ (P37231), THRα (P10827) e THRβ (P10828), restritos a `SINGLE PROTEIN` de *Homo sapiens*.
2. `confidence_score >= 8` (alvo proteico direto).
3. `standard_relation == '='` (descarta `<`, `>`, `~`).
4. Unidade conversível para nM (nM, µM, mM, M, pM).
5. Valor entre 0,001 nM e 10⁷ nM (pActivity de 2 a 12).
6. Descarta registros com `data_validity_comment` problemático (fora da faixa típica, unidade não padrão, possível erro de transcrição ou de autor, duplicata potencial).
7. Exige estrutura química definida (`structure_type == 'MOL'`).
8. Padronização estrutural (`curation.standardize`): limpeza, fragmento parental (remove sais e solventes), normalização, neutralização de carga, tautômero canônico, InChIKey. Rejeita SMILES inválido, elementos não orgânicos, menos de 6 átomos pesados e massa molar fora de 100–1000 Da, com o motivo gravado.
9. Conversão: `pAct = 9 − log10(valor em nM)`. A coluna se chama `pIC50` no dataset de IC50 e `pAct` nas camadas Ki/Kd/EC50.
10. Réplicas: agrupa por InChIKey, usa a **mediana** e descarta o composto se a amplitude entre réplicas passar de 1,0 unidade log.
11. Separa por (alvo, tipo de medida, `assay_type`). Só entram `B` e `F`, com no mínimo 50 registros por grupo.

Cada linha final mantém a proveniência: ID ChEMBL, nome, SMILES original e padronizado, tipo/valor/unidade, pChEMBL, ensaio, alvo, organismo, documento, DOI, PMID e as listas de todos os ensaios e documentos que sustentam o composto.

O `assay_mode` (agonista, antagonista, indefinido) é anotado por expressão regular na descrição do ensaio e fica só como coluna. Nenhum modelo o usa. Por isso os modelos ChEMBL não distinguem agonismo de antagonismo.

BindingDB é usado só como validação externa dos modelos de IC50: ficam os compostos cujo InChIKey não existe no ChEMBL para o mesmo alvo (817 compostos segundo o registro do projeto). O script que calcula as métricas dessa validação não está na pasta.

### 2.3 Tox21 / AC50 (`pipeline_tox21.py`)

A camada Tox21 não é uma regressão de AC50. O que entra no treino é o rótulo ativo/inativo; a pAC50 é calculada e guardada, mas não é usada.

1. Fonte: 12 ensaios qHTS do PubChem fixados na tabela `ASSAYS` (AR e ESR1 em duas linhagens cada, PPARγ e TR em uma, todos em modo agonista e antagonista). Cada ensaio é um dataset.
2. Rótulo: `Active` = 1, `Inactive` = 0, `Inconclusive` = vazio (sai do treino).
3. Pontos censurados: a AC50 existe para todas as substâncias, inclusive inativas. A `pAC50 = 6 − log10(AC50 em µM)` só é preenchida quando o desfecho é `Active`.
4. Citotoxicidade (modos antagonistas): se o contra-ensaio de viabilidade é ativo e sua AC50 é menor ou igual a 3 vezes a AC50 do ensaio-alvo, o "ativo" não é distinguível de morte celular. Ele vira vazio e sai do treino (não vira inativo). A razão de 3 vezes é uma heurística do pipeline, sem calibração formal. Ensaios agonistas não têm contra-ensaio de viabilidade.
5. Autofluorescência: nos ensaios BLA, um "ativo" cujo composto é ativo em qualquer um dos três ensaios de autofluorescência de HEK293 (azul, verde, vermelho) também sai.
6. Réplicas: o Tox21 testa vários SIDs por composto. Resolve por voto majoritário entre as réplicas conclusivas e descarta empates.
7. Estrutura: a mesma padronização do ChEMBL; a chave de composto é o InChIKey.
8. O treino pula ensaios com menos de 30 exemplos na classe minoritária. Por isso `TOX21_TR_agonista_GH3` (27 ativos) não tem modelo e TR só existe em modo antagonista.

### 2.4 Pontos de atenção na camada de dados

1. **Estereoquímica** é registrada (`undefined_stereo` conta centros quirais não definidos), mas não é usada como filtro. Enantiômeros têm InChIKeys diferentes e entram como compostos distintos.
2. **Ki e EC50 compartilham compostos.** No mesmo receptor, 205 dos 271 compostos de `AR_EC50_F` também estão em `AR_Ki_B`; em PPARγ, 177 dos 1.342 compostos de `PPARG_EC50_F` estão em `PPARG_Ki_B`. As medidas são distintas, mas as estruturas se repetem, então os votos dessas camadas não são independentes.
3. **ChEMBL e Tox21 quase não se sobrepõem.** Só 0,4% a 4% dos compostos de cada dataset ChEMBL (Ki/EC50/Kd) aparecem nos ensaios Tox21 do mesmo receptor (máximo 13,8%, `PPARG_Kd_B`). Os dois mundos químicos são diferentes: ligantes de química medicinal contra químicos ambientais.
4. **Divergência a investigar.** Os CSVs em `datasets_tox21_por_ensaio/` não batem exatamente com o `LEIAME_modelos_tox21.csv`. Em `TOX21_AR_agonista_HEK293T` o CSV tem 4.890 compostos e 424 ativos (8,7%), enquanto o LEIAME registra n = 4.844 e fração ativa de 0,075 (cerca de 363 ativos). Não verifiquei a causa; é compatível com uma versão anterior da curadoria.

## 3. Modelos

### 3.1 Features (`features.py`)

2.180 colunas: ECFP4 por contagem (raio 2, 2.048 posições) mais 132 descritores 2D do RDKit (todos exceto os `fr_*`). Descritores com NaN ou infinito viram 0. O ECFP4 binário da mesma molécula é usado só para o domínio de aplicabilidade. O número 2.180 depende da versão do RDKit: com outra versão o preditor falha com `Feature shape mismatch`.

### 3.2 Protocolo de treino (igual nos três scripts de modelagem)

1. Divisão externa 80/20 por scaffold de Bemis-Murcko (`features.scaffold_split`).
2. No conjunto de treino, para cada algoritmo do "zoo", grade pequena de hiperparâmetros com CV 5-fold.
3. Escolha do melhor algoritmo: Q² de CV para regressão, ROC-AUC de CV para classificação (Tox21 também usa ROC-AUC).
4. Avaliação no teste externo por scaffold (R², RMSE e MAE; ROC-AUC, PR-AUC, acurácia balanceada e MCC).
5. y-scrambling no melhor regressor (Q² com rótulos embaralhados deve ficar perto de zero ou negativo).
6. Domínio de aplicabilidade e calibração conformal (3.4).
7. Reajuste do melhor modelo em 100% dos dados curados. Esse é o modelo de produção. As métricas continuam vindo do protocolo com teste externo.

| Tarefa | Algoritmos testados |
|---|---|
| Regressão | Ridge, SVR (RBF), RandomForest (500 árvores), XGBoost (600 árvores) |
| Classificação ChEMBL | Logistic, SVC (RBF), RandomForest (balanceado), XGBoost |
| Classificação Tox21 | Logistic (balanceado), RandomForest (balanceado), XGBoost |

### 3.3 Regressão e classificação do mesmo endpoint

- A regressão prevê `pAct` (contínuo).
- A classificação prevê `pAct >= 6,0` (IC50, Ki, Kd ou EC50 de até 1 µM). O rótulo sai do mesmo valor medido, então **a classificação não usa dados novos**: é o mesmo conjunto de compostos, com o rótulo binarizado.
- Uma classificação só é treinada se a fração de ativos estiver entre 5% e 95% e a classe minoritária tiver pelo menos 10 compostos no treino. `AR_EC50_F` (96% ativos) e `THRB_Ki_B` (96% ativos) não têm classificador; é o motivo de 32 modelos em `models_layers/` e não 34.
- Na camada Tox21 só existe classificação, e o rótulo é a chamada ativo/inativo do depositante, não um corte de potência.

### 3.4 Domínio, conformal e limiar

- **Domínio de aplicabilidade:** o limiar é o percentil 5 da similaridade de Tanimoto ao vizinho mais próximo dentro do próprio treino. Uma molécula nova está no domínio se sua similaridade máxima ao treino for maior ou igual a esse limiar. Valores: 0,29 a 0,63 nos modelos ChEMBL e 0,265 a 0,292 nos de Tox21. Por construção, 5% do próprio treino fica fora, e 0,27 é uma similaridade baixa, então o domínio Tox21 é permissivo.
- **Conformal indutivo (regressão), α = 0,20:** divide o treino em ajuste (75%) e calibração (25%), toma o quantil dos resíduos absolutos e grava `intervalo_inf` e `intervalo_sup`. A cobertura nominal é 80%; a cobertura medida em scaffolds novos fica em `cobertura_scaffold`. O quantil é calibrado num modelo ajustado só com 75% do treino, enquanto o modelo de produção é reajustado com 100%.
- **Limiar de decisão da classificação:** em `models_final`, `conformal.py` calcula o limiar que maximiza o MCC em CV, e ele só foi adotado onde melhorou o teste externo (apenas `PPARG_F`, limiar 0,57; nos demais fica 0,5). Em `models_layers` o limiar é 0,5 fixo. Nos modelos Tox21 o limiar otimizado é sempre adotado, e vai de 0,11 a 0,40 porque as classes são muito desbalanceadas.
- **`nao_implantavel`:** regressão em `models_layers` se n < 200 ou Q² de CV < 0,3; classificação em `models_layers` se n < 200; em `models_final` a marca vem de `LEIAME_calibracao.csv`. Esses modelos são gravados, mas o consenso os ignora.

### 3.5 Como a divisão por scaffold se comporta (verificado)

`scaffold_split` percorre os grupos de scaffold do maior para o menor e manda cada grupo para o teste se ele couber no orçamento de 20%. O docstring diz o contrário ("scaffolds maiores vão para treino"). Rodei a função nos 17 datasets de Ki/Kd/EC50 da pasta:

- O teste fica com as famílias químicas mais populosas e o treino com o restante. Exemplos: `AR_Ki_B` tem 171 compostos de teste vindos de só 4 scaffolds; `PPARG_EC50_B` tem 362 compostos vindos de 13 scaffolds. Em todos os datasets o teste tem entre 2 e 15 scaffolds.
- O tamanho médio do grupo por scaffold vai de 2,7 a 42,8 compostos no teste (mediana 16,5) e de 1,5 a 4,3 no treino. Entre os 5 maiores scaffolds de cada dataset, de 1 a 5 estão no teste.

Consequência: o teste externo é um teste de extrapolação para as famílias químicas dominantes, dependente de poucos scaffolds e pessimista. Isso é compatível com a queda de R² entre CV e teste (a CV é aleatória dentro do treino) e com a cobertura conformal abaixo do nominal. As métricas de teste devem ser lidas como limite inferior, com alta variância. O modelo de produção, reajustado em 100%, contém essas famílias. Não alterei nada; vale decidir qual era a intenção e, se necessário, trocar a regra (grupos grandes para o treino) e reavaliar. Verifiquei só os datasets de Ki/Kd/EC50; os de IC50 usam a mesma função, mas seus CSVs não estão na pasta.

### 3.6 Modelos e desempenho gravados nos LEIAME

**ChEMBL IC50 (`models_final/`, 9 endpoints, 18 modelos).** Cobertura conformal em scaffolds novos; MCC e acurácia balanceada no teste por scaffold com limiar 0,5.

| endpoint | n | reg. | clf. | cobertura conformal (scaffold) | MCC teste | BalAcc teste | limiar | não implantável |
|---|---|---|---|---|---|---|---|---|
| AR_B | 1807 | XGBoost | XGBoost | 0,67 | 0,39 | 0,69 | 0,50 | não |
| AR_F | 618 | RandomForest | RandomForest | 0,79 | 0,00 | 0,50 | 0,50 | sim |
| ESR1_B | 2700 | XGBoost | XGBoost | 0,77 | 0,48 | 0,75 | 0,50 | não |
| ESR1_F | 156 | RandomForest | RandomForest | 0,65 | 0,83 | 0,88 | 0,50 | sim |
| ESR2_B | 1548 | XGBoost | RandomForest | 0,56 | 0,39 | 0,75 | 0,50 | não |
| PPARG_B | 1437 | XGBoost | Logistic | 0,72 | 0,71 | 0,81 | 0,50 | não |
| PPARG_F | 228 | Ridge | SVC | 0,74 | 0,82 | 0,91 | 0,57 | não |
| THRA_B | 220 | RandomForest | RandomForest | 0,70 | 0,28 | 0,58 | 0,50 | sim |
| THRB_B | 396 | XGBoost | XGBoost | 0,89 | 0,49 | 0,78 | 0,50 | não |

**ChEMBL Ki, Kd e EC50 (`models_layers/`, 17 endpoints, 32 modelos).** Q² de CV do melhor regressor e Q² médio com rótulos embaralhados.

| endpoint | n | n_scaffolds | frac. ativos | regressor | Q² CV | Q² scrambling | cob. conformal | limiar AD | classificador | não implantável |
|---|---|---|---|---|---|---|---|---|---|---|
| AR_EC50_B | 232 | 78 | 0,77 | XGBoost | 0,65 | -0,36 | 0,78 | 0,40 | XGBoost | não |
| AR_EC50_F | 271 | 52 | 0,96 | XGBoost | 0,28 | -0,43 | 0,50 | 0,53 | sem classificador | sim |
| AR_Ki_B | 855 | 238 | 0,84 | RandomForest | 0,65 | -0,17 | 0,68 | 0,58 | RandomForest | não |
| ESR1_EC50_B | 661 | 277 | 0,75 | RandomForest | 0,66 | -0,10 | 0,64 | 0,29 | RandomForest | não |
| ESR1_EC50_F | 241 | 71 | 0,73 | Ridge | 0,47 | -0,58 | 0,81 | 0,42 | Logistic | não |
| ESR1_Ki_B | 469 | 188 | 0,82 | XGBoost | 0,66 | -0,35 | 0,57 | 0,48 | SVC | não |
| ESR2_EC50_B | 290 | 98 | 0,74 | XGBoost | 0,40 | -0,37 | 0,67 | 0,41 | RandomForest | não |
| ESR2_EC50_F | 254 | 76 | 0,75 | XGBoost | 0,66 | -0,35 | 0,65 | 0,44 | RandomForest | não |
| ESR2_Ki_B | 351 | 135 | 0,76 | RandomForest | 0,69 | -0,18 | 0,77 | 0,49 | RandomForest | não |
| PPARG_EC50_B | 1808 | 672 | 0,49 | XGBoost | 0,65 | -0,24 | 0,68 | 0,58 | XGBoost | não |
| PPARG_EC50_F | 1342 | 576 | 0,70 | XGBoost | 0,55 | -0,30 | 0,70 | 0,63 | RandomForest | não |
| PPARG_Kd_B | 65 | 36 | 0,52 | Ridge | 0,73 | -0,50 | 1,00 | 0,30 | Logistic | sim |
| PPARG_Ki_B | 466 | 218 | 0,68 | XGBoost | 0,61 | -0,36 | 0,66 | 0,50 | RandomForest | não |
| THRA_EC50_B | 111 | 41 | 0,72 | RandomForest | 0,47 | -0,14 | 0,64 | 0,62 | SVC | sim |
| THRA_Ki_B | 80 | 28 | 0,82 | SVR | 0,36 | -0,34 | 0,81 | 0,54 | XGBoost | sim |
| THRB_EC50_B | 134 | 46 | 0,86 | XGBoost | 0,68 | -0,33 | 0,81 | 0,62 | RandomForest | sim |
| THRB_Ki_B | 104 | 44 | 0,96 | XGBoost | 0,54 | -0,40 | 0,57 | 0,57 | sem classificador | sim |

**Tox21 (`models_tox21/`, 11 modelos de classificação).**

| endpoint | n | frac. ativos | n_scaffolds | modelo | limiar | limiar AD |
|---|---|---|---|---|---|---|
| TOX21_AR_agonista_HEK293T | 4844 | 0,075 | 1462 | RandomForest | 0,40 | 0,278 |
| TOX21_AR_agonista_MDA | 6579 | 0,043 | 2008 | RandomForest | 0,36 | 0,289 |
| TOX21_AR_antagonista_HEK293T | 4838 | 0,147 | 1600 | XGBoost | 0,24 | 0,288 |
| TOX21_AR_antagonista_MDA | 4625 | 0,084 | 1502 | XGBoost | 0,15 | 0,274 |
| TOX21_ESR1_agonista_HEK293T | 6964 | 0,044 | 2182 | RandomForest | 0,26 | 0,292 |
| TOX21_ESR1_agonista_VM7 | 5157 | 0,148 | 1524 | RandomForest | 0,35 | 0,286 |
| TOX21_ESR1_antagonista_HEK293T | 4870 | 0,148 | 1642 | XGBoost | 0,11 | 0,273 |
| TOX21_ESR1_antagonista_VM7 | 4816 | 0,092 | 1601 | RandomForest | 0,20 | 0,273 |
| TOX21_PPARG_agonista_HEK293 | 5776 | 0,030 | 1750 | RandomForest | 0,40 | 0,286 |
| TOX21_PPARG_antagonista_HEK293 | 4522 | 0,100 | 1468 | XGBoost | 0,13 | 0,277 |
| TOX21_TR_antagonista_GH3 | 3457 | 0,222 | 1237 | XGBoost | 0,28 | 0,265 |

## 4. Consenso (`predict_edc2.py`)

### 4.1 Do modelo ao voto

Cada um dos 61 modelos vira uma linha por molécula em `_detalhado.csv`, com um voto binário `ativo`:

| Modelo | Voto "ativo" |
|---|---|
| Regressão | `valor_previsto >= 6,0` (o intervalo conformal é gravado, mas não entra no voto) |
| Classificação | `prob_ativo >= limiar` do modelo |

Antes da contagem, cada linha recebe `no_dominio` (similaridade máxima ao treino do próprio modelo maior ou igual ao limiar de domínio) e `nao_implantavel`. Só contam as linhas com `conta = no_dominio & ~nao_implantavel`.

### 4.2 Do voto ao alerta

Para cada par molécula-receptor, sobre as linhas que contam:

- `camadas_avaliadas` = número de modelos distintos (regressão e classificação do mesmo endpoint contam separadamente);
- `camadas_ativas` = votos ativos; `evidencia_farmacologica` = votos ativos de modelos ChEMBL; `evidencia_tox21` = votos ativos de modelos Tox21.

| Alerta | Critério |
|---|---|
| alto | `camadas_ativas >= 2`, com `evidencia_farmacologica >= 1` e `evidencia_tox21 >= 1` |
| moderado | `fracao_ativa >= 0,5` e `camadas_avaliadas >= 3` |
| baixo | `camadas_ativas >= 1` |
| sem alerta | nenhuma camada ativa |

As regras são avaliadas nessa ordem, e vale a primeira que se aplica.

### 4.3 Limites estruturais do alerta por receptor

O número de camadas possíveis depende de quantos modelos implantáveis existem por receptor. Os modelos Tox21 de TR entram sob o rótulo `TR`, separado de THRA e THRB, e não existe ensaio Tox21 de ERβ. Calculado a partir dos 61 modelos carregados:

| receptor (no resumo) | modelos ChEMBL implantáveis | modelos Tox21 | camadas possíveis | alerta máximo atingível |
|---|---|---|---|---|
| AR | 6 | 4 | 10 | alto |
| ESR1 | 8 | 4 | 12 | alto |
| ESR2 | 8 | 0 | 8 | moderado |
| PPARG | 10 | 2 | 12 | alto |
| THRA | 0 | 0 | 0 | nenhum (não aparece no resumo) |
| THRB | 2 | 0 | 2 | baixo |
| TR | 0 | 1 | 1 | baixo |

Consequências: "alto" só é atingível para AR, ESR1 e PPARG. ESR2 chega no máximo a "moderado". THRβ chega no máximo a "baixo" (dois modelos implantáveis, ambos do mesmo endpoint de IC50). TR, com um único modelo, chega no máximo a "baixo". THRα não tem nenhum modelo ChEMBL implantável e não aparece no resumo. Nas duas execuções em `resultados/` os níveis observados respeitam esses tetos. Portanto "sem alerta" ou "baixo" em ERβ, THR e TR não tem o mesmo peso que em AR, ESR1 e PPARG.

### 4.4 Pontos de atenção do consenso

1. **Dupla contagem do mesmo endpoint.** Regressão e classificação de um endpoint usam os mesmos compostos e as mesmas features, mas contam como duas camadas. Nas execuções `demo` e `teste` são 80 e 185 camadas avaliadas contra 67 e 147 endpoints distintos. Impacto medido nessas duas execuções: colapsar para um voto por endpoint (ativo só se ambos concordam) altera 1 de 22 e 0 de 46 alertas, então o efeito é pequeno, mas o mecanismo existe (por exemplo, Tamoxifen em ESR2 passa de "moderado" para "baixo").
2. **"Ativo" tem definições diferentes por camada.** Nas camadas ChEMBL significa potência de até 1 µM. No Tox21 significa curva ativa na faixa testada, com pAC50 mediana entre 4,4 e 5,9 nos 11 ensaios modelados. Nos CSVs da pasta (ver a divergência em 2.4), só cerca de 15% dos ativos Tox21 teriam pAC50 maior ou igual a 6 (de 7% a 48% por ensaio). Além disso, os limiares de probabilidade do Tox21 (0,11 a 0,40) são bem mais permissivos que o 0,5 do ChEMBL. Um voto Tox21 ativo e um voto ChEMBL ativo não são equivalentes.
3. **Sem direção do efeito.** O alerta soma votos agonistas e antagonistas do Tox21 e votos ChEMBL sem modo. Ele significa "há evidência de interação com o receptor", não agonismo ou antagonismo.
4. **Votos de peso igual, com desempenhos diferentes.** Fora do treino, a classificação generaliza melhor que a regressão (ROC-AUC de 0,76 a 0,94 em 8 de 9 endpoints de IC50, contra R² de teste de −0,26 a 0,66). O voto da regressão pesa o mesmo que o da classificação.
5. **Domínio permissivo no Tox21** (similaridade mínima de cerca de 0,27) e `cobertura_scaffold` abaixo do nominal (leitura em 3.5). Isso reduz a confiança nos votos de moléculas distantes do treino.

### 4.5 Como ler um resultado

1. Veja `alerta` em `_resumo.csv` junto com `camadas_avaliadas` e o teto do receptor (seção 4.3).
2. Abra `_detalhado.csv`, filtre `conta = True` e confira quais modelos votaram. Use `detalhe` para ver o modo Tox21 (agonista ou antagonista) e `tanimoto_viz` para ver a distância ao treino.
3. Trate `valor_previsto` e o intervalo conformal como ordenação e ordem de grandeza, e `prob_ativo` como o indicador principal de atividade.
4. Ausência de alerta não é evidência de ausência de efeito: mecanismos fora dos receptores modelados (aromatase, esteroidogênese) e receptores sem cobertura (IR) não são detectados.

## 5. Decisões em aberto

1. Dupla contagem: manter, usar um voto por endpoint ou usar só a classificação como voto.
2. Regra do `scaffold_split`: manter a atual ou adotar grupos grandes no treino, e reavaliar as métricas.
3. Direção do efeito: documentar no README, acrescentar coluna de modos Tox21 ativos ao `_resumo.csv` ou modelar por ensaio e modo (seção 9 do relatório).
4. Divergência entre os CSVs do Tox21 e o LEIAME (item 4 da seção 2.4).
