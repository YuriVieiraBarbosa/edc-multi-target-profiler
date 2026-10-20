# Verificação da auditoria de código

Li os trechos citados no `pipeline_edc_qsar` e rodei diagnósticos nos 38 endpoints com dataset na pasta
(9 de IC50, 17 de Ki/Kd/EC50 e 12 do Tox21). Nada na sua pasta foi alterado. O que segue separa o que
confirmei, o que precisa de ressalva e o que a auditoria subestimou.

## 1. Itens que confirmei no código

| Item | Evidência |
|---|---|
| A1 | `extract_bindingdb.py:33` grava `cache/bindingdb_ic50_raw.parquet`; `pipeline_curate.py:85` lê `.pkl`. |
| A2 | `bindingdb_std_lookup.pkl` aparece só em `pipeline_curate.py:86`, numa leitura. Nenhum script grava. |
| A3 | `LEIAME_calibracao.csv` tem 22 colunas; `conformal.py` grava `results/conformal_e_limiares.csv`. As 5 colunas (`limiar_producao`, `ganho_limiar`, `conformal_modo`, `cobertura_scaffold`, `nao_implantavel`) vêm de um script que não está na pasta. |
| A5 | `data/` e `results/` não existem. `finalize_models.py:12,15` e `conformal.py` dependem dos dois. |
| A7 | `rodar_edc.py:28,48` chamam `subprocess` sem `cwd`. |
| B1 | `predict_edc2.py:32` compara o nome da pasta como texto. Já documentado em três lugares. |
| B7 | O `agg` usa `ev.loc[s.index, "ativo"]`, dependente de alinhamento de índice. |
| C1 | O docstring diz "scaffolds maiores vão para treino"; o laço preenche o **teste** primeiro. |
| C2–C5 | `conformal.py:50-57` calibra num modelo ajustado em 75% do treino; `finalize_models.py:26` ajusta em 100%. O limiar de domínio é recalculado sobre 100% (`finalize_models.py:27-28`), e `Q2_cv`/`R2_teste_scaffold` do modelo de 80% são gravados dentro do artefato de produção. |

**C2 fecha uma lacuna minha.** No relatório do teste justo eu medi que a subcobertura conformal
persiste nas duas regras de split, e disse que não sabia a causa. A causa provável está aqui: o quantil
é calibrado num modelo e aplicado a outro, ajustado em mais dados. Isso não foi testado, mas é uma
explicação que o split não explicava.

## 2. Item que precisa de ressalva

**A4 não é bug de código.** `pipeline_layers.py:106` grava em `data_layers/` e `modeling_layers.py:95`
lê de `data_layers/`; o mesmo vale para `data_tox21/`. O código é internamente consistente. A
divergência é que na pasta entregue essas saídas estão como `datasets_camadas_por_endpoint/` e
`datasets_tox21_por_ensaio/`, ou seja, foram renomeadas depois de geradas. A correção é de empacotamento
(um nome só, em um lugar só), não de lógica.

## 3. B4: o fallback existe, mas nunca dispara

`features.py:73-75` cai para divisão aleatória quando o teste guloso fica abaixo de metade do alvo, sem
registrar nada. Confirmei o código. Rodei o laço guloso nos 38 endpoints:

- **o fallback não dispara em nenhum deles.** Em todos, o preenchimento guloso atinge exatamente o
  tamanho-alvo do teste (razão 1,00 entre obtido e alvo, no pior caso).

Então nenhum endpoint foi validado aleatoriamente, e o teste justo do `scaffold_split` que rodei não foi
afetado. Isso rebaixa o B4 de "endpoints podem ter sido validados aleatoriamente" para risco latente:
vale registrar a flag, mas não há número a revisar.

## 4. B5 é pior do que a auditoria descreve

Moléculas sem scaffold de Murcko (acíclicas, ou falha da função) caem todas no grupo `""`. Medi o
tamanho desse pseudo-grupo e onde ele cai no split:

- **No Tox21 ele é enorme:** 18% a 23% de cada ensaio, de 654 a 1.392 compostos. Em 7 dos 12 ensaios é
  o maior grupo do dataset.
- **O efeito é bimodal, e depende só de o grupo caber ou não no orçamento do teste.** Em 15 dos 18
  endpoints com acíclicos, o teste fica com **zero** moléculas acíclicas. Nos outros três, acontece o
  contrário:

| Endpoint | Acíclicos | n do teste | Fração do teste que é acíclica |
|---|---|---|---|
| `TOX21_AR_antagonista_HEK293T` | 944 | 992 | 0,95 |
| `TOX21_TR_antagonista_GH3` | 654 | 712 | 0,92 |
| `PPARG_Kd_B` | 4 | 13 | 0,31 |

Ou seja, nesses dois ensaios do Tox21 o teste externo é, na prática, "todas as moléculas acíclicas", e
nos outros nove é "nenhuma molécula acíclica". As duas situações são ruins por motivos opostos, e a
diferença entre elas é só o acaso de `0,2 × n` ser maior ou menor que o pseudo-grupo.

Isso afeta diretamente o ROC-AUC de teste desses dois endpoints, que entraram nas minhas comparações. A
correção mínima é tratar cada molécula acíclica como seu próprio grupo, e isso muda números.

## 5. Itens do `consenso_v2.py`, que são meus — corrigidos

- **E9** (ordem dependia da entrada): `groupby` passou a `sort=True` e os pares saem ordenados.
- **E10** (colunas ausentes): `carregar()` agora cria `valor_previsto`, `intervalo_inf`, `intervalo_sup`,
  `prob_ativo`, `limiar`, `tanimoto_viz` e `detalhe` quando faltam, então um `_detalhado.csv` sem
  calibração conformal não quebra mais.
- **E11** (O(n²)): um `groupby` por `(nome, receptor, eixo)` substituiu o refiltro por par.
- **E12** (`.iloc[0]` silencioso): agora levanta `ValueError` nomeando o endpoint e apontando a causa
  habitual, que é nome repetido na entrada (o B2). Isso transforma o caminho de corrupção silenciosa do
  B2 em erro explícito, pelo menos na reagregação.

Verifiquei que o conteúdo não mudou: as 126 linhas do resumo e as 666 do arquivo de endpoints têm os
mesmos valores de antes, só em ordem determinística. O `ValueError` foi testado duplicando um nome.

**D3 também é meu em parte.** A regra de alerta está copiada no `ablacao_camadas_painel.py`. Ela existe
ali para reproduzir o `_resumo.csv` e provar que a ablação parte do mesmo cálculo, e essa reprodução é
conferida par a par. Ainda assim, é uma terceira cópia.

## 6. O que eu priorizaria

1. **A3.** Enquanto o script de pós-processamento não for reconstruído, os 18 modelos de IC50 não podem
   ser regerados por nada na pasta. É o único item que impede refazer o que está implantado.
2. **B1 e B2.** São os dois caminhos de resultado errado sem aviso. O B2 agora falha de forma explícita
   na reagregação, mas continua fundindo linhas no `predict_edc2.py`.
3. **B5.** Não está na lista de prioridades da auditoria e, pelos números acima, deveria estar: dois
   endpoints do Tox21 têm teste externo praticamente degenerado.
4. **C2 a C4.** São decisões, mas explicam dois achados já medidos (subcobertura conformal e domínio
   implantado mais largo que o validado). Vale congelar com a explicação registrada.

## 7. Ressalvas desta verificação

- Li os trechos citados, não os 24 arquivos inteiros. Itens de arquitetura e de estilo (D e E) não foram
  auditados um a um; confirmei só os que tocam o `consenso_v2.py`.
- Não testei as correções propostas para os itens **D**, porque todas mudam números e dependem da sua
  decisão.
- O diagnóstico do split usa os datasets que estão na pasta mais o artefato de IC50 do projeto. Se os
  datasets forem regerados com outra versão do RDKit, a canonicalização tautomérica muda os InChIKeys e
  o diagnóstico precisa ser refeito (item C9, que não verifiquei).
