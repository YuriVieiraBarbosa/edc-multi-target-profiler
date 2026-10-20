# Teste justo: o efeito de pôr grandes grupos de scaffold no treino

Cópia de trabalho no sandbox; nada em `~/Downloads/pipeline_edc_qsar` foi alterado. O pipeline atual continua como está.

## 1. Resumo

Pergunta: é possível melhorar as métricas adotando os grandes grupos de scaffold no treino? Respondi medindo os três treinos no **mesmo conjunto de teste**, o que a reavaliação anterior não fazia (cada regra de split gerava um teste diferente, e as métricas não eram comparáveis).

1. **Há um ganho pequeno, e ele vem da quantidade de dados, não dos grupos grandes em si.** Com os grandes grupos no treino (condição B, 80% dos dados), a mediana do RMSE cai 0,018 (IC 95%: 0,003 a 0,030; 19 de 26 endpoints melhoram) e o R² sobe 0,021. Em relação à mediana do RMSE de 0,87, isso é cerca de 2%.
2. **A composição sozinha não muda nada.** Quando o treino com grandes grupos é reduzido ao tamanho do treino sem eles (condição C), as diferenças somem: RMSE −0,007 (IC 95%: −0,019 a 0,016), R² −0,003 (IC 95%: −0,018 a 0,012).
3. **Classificação:** o ROC-AUC dos endpoints ChEMBL sobe 0,007 com B (IC 95%: −0,001 a 0,015) e o do Tox21, 0,006 (IC 95%: −0,003 a 0,018). Os dois intervalos incluem zero.
4. **O modelo de produção não muda.** Ele já é treinado com 100% dos dados, grupos grandes inclusos. A regra de split só altera como o desempenho é estimado.
5. **O salto de R² da reavaliação anterior (mediana 0,02 para 0,41) não é ganho de modelo.** No mesmo teste, o ganho de R² é de 0,02. O restante é compatível com a composição do teste antigo, que tinha faixa de atividade restrita (medida na reavaliação anterior); não refiz essa decomposição aqui.

## 2. Desenho

- **Grupos grandes (L):** os grupos de scaffold que `features.scaffold_split` manda para o teste (gulosos, do maior ao menor, até 20% dos compostos). **S** é o restante.
- **Teste comum (T):** cerca de 20% dos compostos, sorteados por grupo de scaffold dentro de S, com semente. Nenhum composto de T aparece em nenhum treino. Mediana por grupo de dados: IC50: 124 compostos de 50 scaffolds; Ki/Kd/EC50: 54 compostos de 21 scaffolds; Tox21: 991 compostos de 750 scaffolds. **R** = S sem T.
- **Condições de treino**, todas avaliadas em T:
  - **A:** treino = R (grandes grupos fora; 60% dos dados).
  - **B:** treino = R ∪ L (grandes grupos dentro; 80% dos dados).
  - **C:** subamostra de R ∪ L com o tamanho de R (isola a composição da quantidade de dados).
- **Modelos:** o algoritmo e os hiperparâmetros que o protocolo escolheu no regime antigo, fixos nas três condições (não reexecutei a seleção por CV). 26 endpoints ChEMBL (9 de IC50 e 17 de Ki/Kd/EC50; 24 com classificador) e 11 endpoints Tox21 (só classificação).
- **Sorteios do teste:** 4 para o ChEMBL (sementes 11 a 14) e 3 para o Tox21 (sementes 11 a 13). Por endpoint, uso a mediana entre sorteios.
- **Critérios fixados antes de olhar os resultados:** unidade de análise = endpoint; teste de Wilcoxon pareado; IC 95% por bootstrap entre endpoints; exclusão de qualquer par endpoint-sorteio com R² < −5 (falha numérica). Nenhuma linha foi excluída por esse critério.

| Grupo | Condição | n de treino (mediana) | Fração do treino vinda de grupos grandes | Tanimoto máximo do teste ao treino (mediana) |
|---|---|---|---|---|
| IC50 | A | 370 | 0,00 | 0,652 |
| IC50 | B | 494 | 0,25 | 0,659 |
| IC50 | C | 370 | 0,25 | 0,643 |
| Ki/Kd/EC50 | A | 163 | 0,00 | 0,635 |
| Ki/Kd/EC50 | B | 217 | 0,25 | 0,644 |
| Ki/Kd/EC50 | C | 163 | 0,25 | 0,635 |
| Tox21 | A | 2975 | 0,00 | 0,378 |
| Tox21 | B | 3966 | 0,25 | 0,404 |
| Tox21 | C | 2975 | 0,25 | 0,386 |

## 3. Resultados

Δ = condição − referência. Para o RMSE, Δ negativo é melhor; para R² e ROC-AUC, positivo é melhor. As colunas de contagem mostram quantos endpoints têm diferença positiva e negativa.

| Métrica | Comparação | Endpoints | Δ mediano | IC 95% | Diferenças > 0 / < 0 | p (Wilcoxon) |
|---|---|---|---|---|---|---|
| RMSE (regressão) | B - A | 26 | −0,018 | −0,030 a −0,003 | 7 / 19 | 0,007 |
| RMSE (regressão) | C - A | 26 | −0,007 | −0,019 a 0,016 | 12 / 14 | 0,727 |
| RMSE (regressão) | B - C | 26 | −0,019 | −0,038 a −0,001 | 8 / 18 | 0,056 |
| R² (regressão) | B - A | 26 | 0,021 | 0,002 a 0,036 | 20 / 6 | 0,024 |
| R² (regressão) | C - A | 26 | −0,003 | −0,018 a 0,012 | 12 / 14 | 0,784 |
| R² (regressão) | B - C | 26 | 0,027 | 0,016 a 0,039 | 20 / 6 | 0,038 |
| ROC-AUC (ChEMBL) | B - A | 24 | 0,007 | −0,001 a 0,015 | 16 / 8 | 0,031 |
| ROC-AUC (ChEMBL) | C - A | 24 | −0,005 | −0,012 a 0,004 | 10 / 14 | 0,422 |
| ROC-AUC (ChEMBL) | B - C | 24 | 0,012 | 0,008 a 0,022 | 18 / 6 | 0,034 |
| ROC-AUC (Tox21) | B - A | 11 | 0,006 | −0,003 a 0,018 | 6 / 5 | 0,123 |
| ROC-AUC (Tox21) | C - A | 11 | −0,010 | −0,018 a 0,002 | 3 / 8 | 0,365 |
| ROC-AUC (Tox21) | B - C | 11 | 0,008 | −0,002 a 0,020 | 8 / 3 | 0,019 |

- **Consistência entre sorteios:** a condição B tem RMSE menor que a A em 66% dos 104 pares endpoint-sorteio, e a C em 45%. No ROC-AUC, B é melhor que A em 62% dos pares do ChEMBL e em 70% dos do Tox21 (C: 45% e 39%).
- **Tamanho do endpoint:** o ganho de B sobre A não se correlaciona de modo detectável com o tamanho do treino (Spearman −0,16 para o RMSE e −0,22 para o ROC-AUC, p > 0,3; n = 26 e 24).
- **O teste fica mais perto do treino só um pouco:** o Tanimoto máximo mediano do teste ao treino passa, por exemplo, de 0,652 para 0,659 em IC50 e de 0,378 para 0,404 em Tox21.

![Ganho por endpoint, condições B e C em relação a A]({{artifact:art_ea568f8f-7232-4974-8e25-fbe9fe645b42}})

## 4. Interpretação

1. Adicionar os grupos grandes ao treino melhora o desempenho em compostos de scaffolds novos, mas pouco: cerca de 0,02 unidade de pAct no RMSE, 0,02 no R² e 0,007 no ROC-AUC dos endpoints ChEMBL.
2. O ganho é explicado pelo aumento de 60% para 80% dos dados: B − C (só quantidade) tem Δ de −0,019 no RMSE, 0,027 no R² e 0,012 no ROC-AUC, e C − A (só composição) fica próximo de zero. Não há evidência de que os grupos grandes sejam mais informativos por composto que os pequenos, embora o IC do contraste C − A ainda admita um efeito de até cerca de 0,02 no RMSE.
3. Para a decisão sobre a regra do `scaffold_split`: adotar grupos grandes no treino é defensável como protocolo de avaliação (treina com 80% e testa em scaffolds menos populosos, parecido com a priorização de compostos novos), mas não é uma alavanca de melhoria do modelo. Ganhos maiores exigiriam outros caminhos, que não testei (por exemplo, seleção de modelo por CV agrupada por scaffold e transformação do descritor `Ipc`).
4. Ao comparar métricas entre regras de split, compare em um teste comum. Os números de R² da regra corrigida e da regra atual não são diretamente comparáveis.

## 5. Limitações

1. **Teste só com scaffolds pequenos.** T é sorteado dentro de S, então a conclusão vale para scaffolds novos e pouco populosos, e não para as famílias químicas dominantes (L). O teste da regra atual (extrapolação para L) não é comparável.
2. **A seleção de algoritmo e de hiperparâmetros usou CV em S, que contém T.** Isso pode inflar levemente as métricas absolutas das três condições. As diferenças entre condições são menos afetadas, porque a configuração é a mesma, mas o efeito não foi quantificado. A configuração não foi reotimizada para a condição B.
3. **Tamanhos diferentes dos usados antes.** O treino da condição A tem 60% dos dados, então as métricas absolutas não são comparáveis com as do relatório anterior (80%).
4. **Poucos sorteios** (4 no ChEMBL, 3 no Tox21), e a quinta semente planejada para o ChEMBL não rodou por um erro de lista de argumentos na fila.
5. **Endpoints não independentes.** IC50, Ki, Kd e EC50 do mesmo receptor compartilham compostos, então os ICs entre endpoints tendem a ser otimistas.
6. **Sem correção para múltiplas comparações** (24 testes no total: 8 métricas e 3 contrastes). Os p isolados são orientativos. Só os contrastes B − A e B − C apontam na mesma direção em várias métricas.
7. **Tox21:** os CSVs da pasta não batem com os modelos implantados (divergência já registrada); este teste vale para os CSVs.
8. **Falhas do Ridge com `Ipc`:** nenhuma falha numérica apareceu aqui, embora 3 regressões usem Ridge no regime antigo. Não investiguei por que o problema, visto antes com a regra corrigida, não se repete neste teste; isso não elimina o risco.

## 6. Arquivos

`fair_test.py` executa o teste por endpoint e `analise_teste_justo.py` gera as tabelas; `fair_longo.csv` traz uma linha por endpoint, condição, tarefa e sorteio; `fair_resumo_pareado.csv` traz todas as comparações, incluindo MAE, Spearman e PR-AUC.

### Todas as comparações

| Métrica | Comparação | Endpoints | Δ mediano | IC 95% | Diferenças > 0 / < 0 | p (Wilcoxon) |
|---|---|---|---|---|---|---|
| RMSE (regressão) | B - A | 26 | −0,018 | −0,030 a −0,003 | 7 / 19 | 0,007 |
| RMSE (regressão) | C - A | 26 | −0,007 | −0,019 a 0,016 | 12 / 14 | 0,727 |
| RMSE (regressão) | B - C | 26 | −0,019 | −0,038 a −0,001 | 8 / 18 | 0,056 |
| MAE (regressão) | B - A | 26 | −0,012 | −0,023 a 0,007 | 10 / 16 | 0,075 |
| MAE (regressão) | C - A | 26 | −0,000 | −0,012 a 0,022 | 13 / 13 | 0,764 |
| MAE (regressão) | B - C | 26 | −0,016 | −0,031 a 0,016 | 10 / 16 | 0,129 |
| R² (regressão) | B - A | 26 | 0,021 | 0,002 a 0,036 | 20 / 6 | 0,024 |
| R² (regressão) | C - A | 26 | −0,003 | −0,018 a 0,012 | 12 / 14 | 0,784 |
| R² (regressão) | B - C | 26 | 0,027 | 0,016 a 0,039 | 20 / 6 | 0,038 |
| Spearman (regressão) | B - A | 26 | 0,018 | 0,000 a 0,026 | 18 / 8 | 0,080 |
| Spearman (regressão) | C - A | 26 | −0,008 | −0,018 a 0,003 | 9 / 17 | 0,280 |
| Spearman (regressão) | B - C | 26 | 0,019 | 0,010 a 0,029 | 20 / 6 | 0,025 |
| ROC-AUC (ChEMBL) | B - A | 24 | 0,007 | −0,001 a 0,015 | 16 / 8 | 0,031 |
| ROC-AUC (ChEMBL) | C - A | 24 | −0,005 | −0,012 a 0,004 | 10 / 14 | 0,422 |
| ROC-AUC (ChEMBL) | B - C | 24 | 0,012 | 0,008 a 0,022 | 18 / 6 | 0,034 |
| ROC-AUC (Tox21) | B - A | 11 | 0,006 | −0,003 a 0,018 | 6 / 5 | 0,123 |
| ROC-AUC (Tox21) | C - A | 11 | −0,010 | −0,018 a 0,002 | 3 / 8 | 0,365 |
| ROC-AUC (Tox21) | B - C | 11 | 0,008 | −0,002 a 0,020 | 8 / 3 | 0,019 |
| PR-AUC (ChEMBL) | B - A | 24 | 0,007 | 0,003 a 0,014 | 20 / 4 | 0,001 |
| PR-AUC (ChEMBL) | C - A | 24 | −0,000 | −0,006 a 0,009 | 12 / 12 | 0,546 |
| PR-AUC (ChEMBL) | B - C | 24 | 0,007 | 0,003 a 0,011 | 19 / 5 | 0,016 |
| PR-AUC (Tox21) | B - A | 11 | 0,012 | −0,004 a 0,030 | 8 / 3 | 0,032 |
| PR-AUC (Tox21) | C - A | 11 | 0,002 | −0,025 a 0,023 | 6 / 5 | 0,966 |
| PR-AUC (Tox21) | B - C | 11 | 0,011 | 0,006 a 0,024 | 10 / 1 | 0,002 |
