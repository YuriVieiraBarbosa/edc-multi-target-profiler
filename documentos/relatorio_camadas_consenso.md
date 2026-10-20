# Camadas paralelas de bioatividade e predição por consenso

Extensão do pipeline QSAR de desregulação endócrina para superar as três
limitações identificadas na primeira entrega: recorte restrito a IC50, ausência
de química ambiental e domínio de aplicabilidade binário.

---

## 1. Decisões tomadas

| Limitação | Decisão | Resultado medido |
|---|---|---|
| Só IC50 — genisteína e rosiglitazona ficavam de fora | Extrair Ki, Kd e EC50 do ChEMBL 37 como **camadas independentes**, nunca fundidas com o IC50 | 17 datasets novos, 7.733 pares composto–endpoint |
| Contaminantes ambientais fora do domínio | Incorporar o **Tox21 qHTS** completo (biblioteca de ~10.500 substâncias) via arquivo FTP do PubChem | 12 datasets, 7226 compostos únicos |
| "Atividade" no qHTS confunde antagonismo com morte celular | Usar os **contra-ensaios de viabilidade pareados** e os controles de autofluorescência | 15% a 55% dos ativos de antagonismo reclassificados |
| Domínio de aplicabilidade binário, sem incerteza por composto | **Predição conformal indutiva** com cobertura nominal declarada e verificada | cobertura empírica 0,57–0,89 em scaffolds novos |
| Classificadores desbalanceados travados em 0,5 | Otimizar o limiar por validação cruzada e **manter 0,5 quando a otimização não transfere** | só PPARγ funcional melhorou (MCC 0,82 → 0,91) |
| Predição de alvo único, sem integração | **Preditor de consenso** com evidência por camada e escore de alerta | sensibilidade 12/15, especificidade 3/3 no painel de referência |

![Expansão em camadas]({{artifact:art_521608b2-0ad8-4598-8310-682fd8724f91}})

---

## 2. Camadas de dados

Cada tipo de medida é um endpoint próprio. A comparabilidade de ensaio continua
sendo o critério que governa tudo: alvo humano SINGLE PROTEIN, `confidence_score ≥ 8`,
relação `=`, unidade conversível a nM, mediana por InChIKey com descarte de
réplicas com amplitude acima de 1 unidade log, e separação entre ensaio de
ligação e funcional.

### 2.1 ChEMBL 37 — Ki, Kd e EC50

15.007 registros brutos, 10.273 após os filtros.

| receptor   | medida   | ensaio    |   registros |   compostos_finais |   pAct_mediana |   descartados_discordantes |
|:-----------|:---------|:----------|------------:|-------------------:|---------------:|---------------------------:|
| PPARγ      | EC50     | ligação   |        2356 |               1808 |           5.96 |                         35 |
| PPARγ      | EC50     | funcional |        1764 |               1342 |           6.52 |                         51 |
| AR         | Ki       | ligação   |        1155 |                855 |           7.21 |                         36 |
| ERα        | EC50     | ligação   |         822 |                661 |           7.28 |                         35 |
| ERα        | Ki       | ligação   |         648 |                469 |           6.99 |                         22 |
| PPARγ      | Ki       | ligação   |         566 |                466 |           6.6  |                         12 |
| ERβ        | Ki       | ligação   |         405 |                351 |           6.92 |                          4 |
| ERβ        | EC50     | ligação   |         411 |                290 |           6.69 |                         17 |
| AR         | EC50     | funcional |         311 |                271 |           7.82 |                          3 |
| ERβ        | EC50     | funcional |         350 |                254 |           6.98 |                         19 |
| ERα        | EC50     | funcional |         293 |                241 |           6.48 |                          6 |
| AR         | EC50     | ligação   |         278 |                232 |           8.01 |                          6 |
| THRβ       | EC50     | ligação   |         178 |                134 |           6.89 |                          7 |
| THRα       | EC50     | ligação   |         149 |                111 |           6.66 |                          4 |
| THRβ       | Ki       | ligação   |         118 |                104 |           8.47 |                          2 |
| THRα       | Ki       | ligação   |          88 |                 80 |           6.97 |                          1 |
| PPARγ      | Kd       | ligação   |          81 |                 65 |           6.02 |                          1 |

### 2.2 Tox21 — AC50 do qHTS

Doze ensaios de receptor nuclear do arquivo FTP do PubChem (bloco
2283001–2284000), com os contra-ensaios de viabilidade pareados e os controles
de autofluorescência em HEK293. A AC50 vem em µM e é reportada pelo PubChem
inclusive para substâncias inativas — só é usada quando o desfecho é *Active*.

| receptor   | modo        | linhagem    |   n_testados |   n_ativos |   n_citotox_confunde |   compostos_validos |   frac_ativos |
|:-----------|:------------|:------------|-------------:|-----------:|---------------------:|--------------------:|--------------:|
| ERα        | agonista    | HEK293T-BLA |        10496 |        594 |                    0 |                6965 |        0.0441 |
| ERα        | antagonista | HEK293T-BLA |        10496 |       1315 |                  198 |                4871 |        0.1482 |
| ERα        | agonista    | VM7-LUC     |        10496 |       1283 |                    0 |                5158 |        0.1483 |
| ERα        | antagonista | VM7-LUC     |        10496 |       1028 |                  316 |                4817 |        0.0924 |
| AR         | agonista    | HEK293T-BLA |        10496 |        625 |                    0 |                4845 |        0.0747 |
| AR         | antagonista | HEK293T-BLA |        10496 |       1915 |                  761 |                4839 |        0.1469 |
| AR         | agonista    | MDA-kb2-LUC |        10496 |        501 |                    0 |                6580 |        0.0433 |
| AR         | antagonista | MDA-kb2-LUC |        10496 |       1196 |                  569 |                4626 |        0.0841 |
| PPARγ      | agonista    | HEK293-BLA  |        10496 |        296 |                    0 |                5777 |        0.0301 |
| PPARγ      | antagonista | HEK293-BLA  |         9315 |        858 |                  198 |                4522 |        0.1    |
| THR*       | agonista    | GH3-LUC     |        10496 |         35 |                    0 |                6729 |        0.0034 |
| THR*       | antagonista | GH3-LUC     |        10496 |       2651 |                 1459 |                3458 |        0.2224 |

Duas regras de curadoria específicas desta camada:

1. **Confundimento por citotoxicidade.** Um composto ativo cuja AC50 de
   viabilidade é menor ou igual a 3× a AC50 do ensaio-alvo não é distinguível
   de morte celular. Esses registros recebem rótulo ausente — **não** são
   reclassificados como inativos, porque a informação sobre o receptor
   simplesmente não existe naquela medida.
2. **Réplicas.** O Tox21 testa várias amostras (SIDs) do mesmo composto, e elas
   discordam. A resolução é por voto majoritário entre as réplicas conclusivas,
   com empate descartado. Aplicar essa regra melhorou o ROC-AUC no teste por
   scaffold de 0,881 para 0,922 em AR agonista (HEK293T) e de 0,841 para 0,884
   em AR agonista (MDA-kb2) — o ganho é remoção de ruído de rótulo, não ajuste.

### 2.3 As camadas não são intercambiáveis

Nos compostos em comum, a correlação entre camadas vai de r = −0,24 a r = 0,90,
com erro absoluto mediano de 0,4 a 0,9 unidades log mesmo nos pares mais
correlacionados. O par mais discordante — AR EC50 funcional contra AR IC50
funcional — mede agonismo contra antagonismo no mesmo tipo de ensaio, e é
justamente onde fundir os dados produziria o pior dano.

| camada_A     | camada_B     |   n_comum |   r_pearson |   rho_spearman |   vies_A_menos_B |   erro_absoluto_mediano |
|:-------------|:-------------|----------:|------------:|---------------:|-----------------:|------------------------:|
| PPARG|EC50|F | PPARG|IC50|B |       400 |       0.809 |          0.812 |             0.01 |                    0.41 |
| PPARG|EC50|B | PPARG|IC50|B |       277 |       0.619 |          0.635 |            -0.24 |                    0.68 |
| AR|IC50|F    | AR|Ki|B      |       265 |       0.568 |          0.557 |            -0.21 |                    0.45 |
| AR|EC50|F    | AR|Ki|B      |       205 |       0.538 |          0.555 |            -0.08 |                    0.52 |
| PPARG|EC50|F | PPARG|Ki|B   |       177 |       0.735 |          0.736 |            -0.23 |                    0.46 |
| AR|IC50|B    | AR|Ki|B      |       176 |       0.567 |          0.636 |            -0.05 |                    0.63 |
| ESR1|EC50|B  | ESR1|IC50|B  |       148 |       0.61  |          0.582 |             0.51 |                    0.91 |
| ESR2|EC50|F  | ESR2|IC50|B  |       134 |       0.862 |          0.861 |             0.24 |                    0.42 |
| AR|EC50|B    | AR|IC50|B    |       118 |       0.763 |          0.734 |            -0.03 |                    0.77 |
| ESR1|IC50|B  | ESR1|Ki|B    |       104 |       0.789 |          0.771 |            -0.46 |                    0.61 |
| PPARG|EC50|F | PPARG|IC50|F |        91 |       0.759 |          0.758 |            -0.11 |                    0.46 |
| ESR1|EC50|F  | ESR1|IC50|B  |        86 |       0.685 |          0.718 |             0.46 |                    0.62 |

O Tox21 quase não se sobrepõe ao ChEMBL: entre 10 e 26 compostos em comum por
par. Isso confirma que ele cobre espaço químico complementar, não redundante.

---

## 3. Modelos

### 3.1 Camadas Ki / Kd / EC50

Mesmo protocolo dos modelos de IC50: split externo por scaffold de
Bemis-Murcko, seleção por CV 5-fold, y-scrambling, refit de produção em 100%
dos dados.

| receptor   | medida   | ensaio    |    n | melhor_modelo   |   Q2_cv |   y_scrambling_Q2_medio |   cobertura_conformal_scaffold | nao_implantavel   |
|:-----------|:---------|:----------|-----:|:----------------|--------:|------------------------:|-------------------------------:|:------------------|
| AR         | EC50     | ligação   |  232 | XGBoost         |   0.647 |                  -0.365 |                          0.783 | False             |
| AR         | EC50     | funcional |  271 | XGBoost         |   0.276 |                  -0.427 |                          0.5   | True              |
| AR         | Ki       | ligação   |  855 | RandomForest    |   0.649 |                  -0.173 |                          0.684 | False             |
| ERα        | EC50     | ligação   |  661 | RandomForest    |   0.66  |                  -0.104 |                          0.636 | False             |
| ERα        | EC50     | funcional |  241 | Ridge           |   0.474 |                  -0.58  |                          0.812 | False             |
| ERα        | Ki       | ligação   |  469 | XGBoost         |   0.662 |                  -0.348 |                          0.574 | False             |
| ERβ        | EC50     | ligação   |  290 | XGBoost         |   0.402 |                  -0.37  |                          0.672 | False             |
| ERβ        | EC50     | funcional |  254 | XGBoost         |   0.661 |                  -0.347 |                          0.647 | False             |
| ERβ        | Ki       | ligação   |  351 | RandomForest    |   0.688 |                  -0.18  |                          0.771 | False             |
| PPARγ      | EC50     | ligação   | 1808 | XGBoost         |   0.646 |                  -0.244 |                          0.682 | False             |
| PPARγ      | EC50     | funcional | 1342 | XGBoost         |   0.551 |                  -0.299 |                          0.698 | False             |
| PPARγ      | Kd       | ligação   |   65 | Ridge           |   0.727 |                  -0.502 |                          1     | True              |
| PPARγ      | Ki       | ligação   |  466 | XGBoost         |   0.609 |                  -0.363 |                          0.656 | False             |
| THRα       | EC50     | ligação   |  111 | RandomForest    |   0.471 |                  -0.144 |                          0.636 | True              |
| THRα       | Ki       | ligação   |   80 | SVR             |   0.361 |                  -0.338 |                          0.812 | True              |
| THRβ       | EC50     | ligação   |  134 | XGBoost         |   0.678 |                  -0.327 |                          0.815 | True              |
| THRβ       | Ki       | ligação   |  104 | XGBoost         |   0.54  |                  -0.397 |                          0.571 | True              |

### 3.2 Classificadores Tox21

| endpoint                       | modelo       |   ROC_AUC_cv |   ROC_AUC_teste |   PR_AUC_teste |   MCC_teste_limiar |   limiar |
|:-------------------------------|:-------------|-------------:|----------------:|---------------:|-------------------:|---------:|
| TOX21_AR_agonista_HEK293T      | RandomForest |        0.844 |           0.922 |          0.773 |              0.692 |     0.4  |
| TOX21_AR_agonista_MDA          | RandomForest |        0.788 |           0.884 |          0.788 |              0.788 |     0.36 |
| TOX21_AR_antagonista_HEK293T   | XGBoost      |        0.834 |           0.782 |          0.274 |              0.274 |     0.24 |
| TOX21_AR_antagonista_MDA       | XGBoost      |        0.84  |           0.832 |          0.375 |              0.343 |     0.15 |
| TOX21_ESR1_agonista_HEK293T    | RandomForest |        0.882 |           0.759 |          0.217 |              0.117 |     0.26 |
| TOX21_ESR1_agonista_VM7        | RandomForest |        0.732 |           0.814 |          0.661 |              0.452 |     0.35 |
| TOX21_ESR1_antagonista_HEK293T | XGBoost      |        0.856 |           0.805 |          0.426 |              0.361 |     0.11 |
| TOX21_ESR1_antagonista_VM7     | RandomForest |        0.862 |           0.833 |          0.444 |              0.404 |     0.2  |
| TOX21_PPARG_agonista_HEK293    | RandomForest |        0.834 |           0.796 |          0.22  |              0     |     0.4  |
| TOX21_PPARG_antagonista_HEK293 | XGBoost      |        0.836 |           0.829 |          0.465 |              0.366 |     0.13 |
| TOX21_TR_antagonista_GH3       | XGBoost      |        0.8   |           0.72  |          0.235 |              0.253 |     0.28 |

O ensaio TR agonista (GH3) foi descartado: 35 ativos em 10.496 substâncias não
sustentam modelo.

### 3.3 Calibração conformal e ponto de operação

| endpoint   | modelo       |   q_absoluto |   cobertura_scaffold | conformal_modo   |   limiar_producao |   ganho_limiar |
|:-----------|:-------------|-------------:|---------------------:|:-----------------|------------------:|---------------:|
| AR_B       | XGBoost      |        0.801 |                0.673 | fixo             |              0.5  |         -0.041 |
| AR_F       | RandomForest |        0.764 |                0.79  | fixo             |              0.5  |          0     |
| ESR1_B     | XGBoost      |        0.837 |                0.772 | fixo             |              0.5  |         -0.063 |
| ESR1_F     | RandomForest |        1.234 |                0.645 | fixo             |              0.5  |          0     |
| ESR2_B     | XGBoost      |        0.925 |                0.565 | fixo             |              0.5  |         -0.06  |
| PPARG_B    | XGBoost      |        0.82  |                0.718 | fixo             |              0.5  |         -0.035 |
| PPARG_F    | Ridge        |        1.025 |                0.739 | fixo             |              0.57 |          0.092 |
| THRA_B     | RandomForest |        0.592 |                0.705 | normalizado      |              0.5  |         -0.003 |
| THRB_B     | XGBoost      |        0.715 |                0.886 | fixo             |              0.5  |          0     |

**A cobertura conformal fica abaixo do nominal.** Com α = 0,20, a cobertura
empírica no teste por scaffold é de 0,57 a 0,89 (mediana 0,72), não 0,80. Isso
não é defeito de implementação: a garantia conformal vale sob permutabilidade,
e um split por scaffold quebra a permutabilidade de propósito — é essa a razão
de usá-lo. O número honesto para interpretar uma predição em química nova é a
cobertura medida em scaffolds novos, registrada em
`cobertura_scaffold` dentro de cada modelo, não o α nominal.

A otimização de limiar por CV **não transferiu**: dos 9 endpoints de IC50, o MCC
no teste por scaffold piorou em 5 (AR ligação, ERα ligação, ERβ ligação, PPARγ
ligação, THRα ligação), ficou inalterado em 3 (AR funcional, ERα funcional, THRβ
ligação) e melhorou em 1. O limiar otimizado só foi adotado em PPARγ funcional. Três endpoints estão marcados como não implantáveis
no próprio arquivo do modelo: AR funcional, ERα funcional e THRα ligação.

---

## 4. Predição por consenso

`predict_edc2.py` avalia cada molécula em todas as camadas disponíveis e devolve
duas tabelas: uma detalhada, com uma linha por camada (valor previsto, intervalo
conformal, probabilidade, similaridade ao vizinho mais próximo do treino e flag
de domínio), e um resumo por receptor com o escore de alerta.

| Nível | Regra |
|---|---|
| alto | pelo menos 2 camadas ativas, sendo ao menos uma farmacológica (ChEMBL) e uma ambiental (Tox21) |
| moderado | pelo menos 3 camadas avaliadas e maioria ativa |
| baixo | evidência isolada, 1 camada ativa |
| sem alerta | nenhuma camada ativa dentro do domínio |

Só entram no escore as camadas em que a molécula está dentro do domínio de
aplicabilidade **e** cujo modelo não está marcado como não implantável.

### Painel de referência

19 compostos de entrada: 9 desreguladores endócrinos reconhecidos, 6 referências
de mecanismo (agonista de ERα, agonista de PPARγ, antagonista de AR, fitoestrógeno,
micotoxina estrogênica e antimicrobiano), 3 controles negativos e 1 fármaco não-EDC.

**18 foram avaliados.** A glicina é rejeitada na padronização por massa molar
abaixo de 100 Da — o filtro de MW do pipeline, aplicado a todas as fontes. O
preditor agora grava as moléculas rejeitadas em um arquivo `_rejeitadas.csv` e
as reporta na saída, em vez de apenas omiti-las.

**Sensibilidade 12/15, especificidade 3/3** (entre os 18 avaliados; dos 4
controles não-EDC do painel, a glicina não chegou a ser avaliada). As duas lacunas que motivaram esta
fase foram fechadas: genisteína passa a ERα alto (pelas camadas Ki e EC50) e
rosiglitazona a PPARγ alto. DDT, metoxicloro, vinclozolina e triclosan, antes
fora do domínio, passam a ser recuperados pela camada Tox21.

Os três falsos negativos — atrazina, DBP e DEHP — merecem leitura cuidadosa:
os três são **medidos inativos** nos próprios ensaios Tox21 desses cinco
receptores. Atrazina age por indução de aromatase, DEHP pelo metabólito MEHP e
DBP por interferência na esteroidogênese fetal. Não são falhas de modelo; são
mecanismos fora do escopo de receptor modelado. Ampliar a cobertura para eles
exige endpoints de aromatase (CYP19A1) e de esteroidogênese (H295R), não
modelos melhores dos mesmos receptores.

---

## 5. Limitações que permanecem

1. **ERβ não tem camada ambiental.** O Tox21 não inclui ensaio de ERβ. A
   seletividade α/β continua apoiada só em dados farmacológicos.
2. **THR sem resolução de subtipo.** O ensaio Tox21 de TR é em célula GH3 de
   rato e não distingue α de β; nas tabelas ele aparece como `THR*` e nunca é
   somado aos endpoints humanos THRα e THRβ.
3. **PPARγ agonista no Tox21 tem prevalência de 3%.** O MCC no teste por
   scaffold é 0, apesar do ROC-AUC de 0,80: o modelo ordena, mas não classifica.
   Use a probabilidade como ranqueamento, não como decisão.
4. **A cobertura conformal é menor que a nominal em química nova** (seção 3.3).
5. **Cortellis continua fora.** Produto licenciado da Clarivate, sem API
   pública. Uma exportação CSV entra na mesma estrutura de proveniência.
6. **Nada aqui substitui ensaio in vitro** nem constitui avaliação regulatória
   de risco. O escore de alerta é ferramenta de priorização.

---

## 6. Reprodutibilidade

| Arquivo | Conteúdo |
|---|---|
| `pipeline_layers.py` | extração e curadoria das camadas Ki / Kd / EC50 do ChEMBL 37 |
| `pipeline_tox21.py` | curadoria Tox21: citotoxicidade, autofluorescência, réplicas |
| `modeling_layers.py` | modelos das camadas Ki / Kd / EC50 |
| `modeling_tox21.py` | classificadores Tox21 com limiar otimizado |
| `conformal.py` | calibração conformal e ponto de operação dos modelos de IC50 |
| `predict_edc2.py` | preditor de consenso multi-camada |

Uso:

```bash
python predict_edc2.py --input moleculas.csv --smiles-col smiles --name-col nome
```

Ambiente: Python 3.11, RDKit 2026.03.6, scikit-learn 1.9.1, XGBoost. Semente
fixa (42) em todos os splits e modelos.
