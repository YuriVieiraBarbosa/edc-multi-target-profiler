# Consenso v2: relatório de dois eixos

Implementação do passo 1 da arquitetura proposta: **mudanças só de agregação**. Nenhum modelo foi
retreinado ou alterado, nenhum voto individual mudou, e nada em `~/Downloads/pipeline_edc_qsar` foi
tocado. O `consenso_v2.py` lê o `<saida>_detalhado.csv` que o seu `predict_edc2.py` já produz e
reescreve apenas como esses votos são resumidos.

## 1. O que mudou

| | `_resumo.csv` atual (v1) | Consenso v2 |
|---|---|---|
| Unidade de evidência | o modelo (regressão e classificação do mesmo endpoint contam 2 vezes) | o endpoint (1 evidência) |
| Papel de cada tarefa | as duas votam igual | a classificação decide a positividade; a regressão dá a potência com intervalo conformal |
| Agregação | um nível único, somando ChEMBL e Tox21 | dois eixos, nunca somados: farmacológico e triagem |
| Dentro do eixo | fração de votos ativos | melhor evidência e quantos endpoints a corroboram |
| Modelo fora do domínio | desaparece e dilui a fração | estado explícito `nao_avaliado`, com a cobertura ao lado |
| Direção do efeito | ausente | rótulo a partir do modo dos ensaios Tox21 |
| Pares sem cobertura | sem linha no arquivo | linha com `sem cobertura` |

Por que a classificação decide: fora do treino ela generaliza melhor (ROC-AUC de teste de 0,76 a 0,94
nos endpoints de IC50, contra R² de teste de −0,26 a 0,66). A regressão continua no relatório como
ordem de grandeza da potência, com o intervalo conformal, e não como um segundo voto.

Por que a melhor evidência em vez da fração: os endpoints não são independentes. Ki e EC50 compartilham
compostos (205 dos 271 de `AR_EC50_F` estão em `AR_Ki_B`), e os 4 ensaios Tox21 de AR e de ERα são
2 modos em 2 linhagens. A fração dá peso implícito ao receptor que tem mais datasets.

## 2. Estados e leitura

Cada eixo recebe um estado: `positivo`, `negativo`, `nao_avaliado` (há modelo implantável, mas a
molécula está fora do domínio) ou `sem_modelo` (não existe modelo implantável para aquele receptor
naquele eixo). A leitura conjunta nomeia as 9 combinações possíveis, sem somar nada.

## 3. Cobertura por receptor

Este é o teto de cada receptor, igual para qualquer molécula. Ele explica por que um resultado
negativo em ERβ ou em THRβ vale menos que em AR ou ERα.

| Receptor | Endpoints farmacológicos | Ensaios Tox21 | Medidas | Modos Tox21 |
|---|---|---|---|---|
| AR | 3 | 4 | EC50, IC50, Ki | agonista, antagonista |
| ESR1 | 4 | 4 | EC50, IC50, Ki | agonista, antagonista |
| ESR2 | 4 | 0 | EC50, IC50, Ki | — |
| PPARG | 5 | 2 | EC50, IC50, Ki | agonista, antagonista |
| THRB | 1 | 0 | IC50 | — |
| TR | 0 | 1 | — | antagonista |

THRα não aparece: os 6 modelos ChEMBL dele são todos não implantáveis e não há ensaio Tox21 com esse
rótulo. ERβ não tem ensaio Tox21, e o rótulo `TR` do Tox21 não resolve α de β.

## 4. Comparação com o resumo atual, no painel de referência

Rodei o seu `predict_edc2.py` sem alterações no `painel_ref.csv` (18 moléculas) e reagreguei o
`_detalhado.csv`. O v1 produziu 83 linhas, e o v2 produz 126 pares molécula × receptor
(43 deles não tinham linha no v1).

| Leitura v2 | alto | moderado | baixo | sem alerta | (sem linha no v1) |
|---|---|---|---|---|---|
| convergente negativo | 0 | 0 | 0 | 1 | 0 |
| convergente positivo | 10 | 0 | 0 | 0 | 0 |
| divergente: farmacologico positivo, triagem negativo | 0 | 0 | 1 | 0 | 0 |
| divergente: triagem positivo, farmacologico negativo | 0 | 1 | 3 | 0 | 0 |
| negativo só no eixo de triagem (cobertura parcial) | 0 | 0 | 0 | 34 | 0 |
| negativo só no eixo farmacologico (cobertura parcial) | 0 | 0 | 0 | 5 | 0 |
| sem cobertura | 0 | 0 | 0 | 0 | 43 |
| só farmacologico positivo (sem cobertura de triagem) | 0 | 4 | 2 | 0 | 0 |
| só triagem positivo (sem cobertura farmacologica) | 0 | 5 | 17 | 0 | 0 |

O que a tabela mostra:

1. **Os 10 pares "alto" do v1 viram todos "convergente positivo".** O nível mais forte do v1 já era,
   de fato, convergência entre os dois eixos.
2. **Os 10 pares "moderado" se separam em três situações distintas**, que o nível único juntava:
   5 são só triagem positiva sem nenhuma cobertura farmacológica (Metoxicloro, Vinclozolina e
   Triclosan em AR; Triclosan e Flutamida em ERα), 4 são só farmacológico positivo sem cobertura de
   triagem (os quatro de ERβ, que não tem ensaio Tox21), e 1 é divergência (Zearalenona em ERα, com
   triagem positiva e farmacológico negativo).
3. **Dos 40 pares "sem alerta" do v1, só 1 é negativo nos dois eixos.** Os outros 39 são negativos em
   um eixo com o outro sem cobertura. Na prática, "sem alerta" quase nunca significava que as duas
   linhas de evidência disseram não.
4. **O viés de cobertura fica visível em vez de embutido.** Os 5 moderados de AR e ERα sem cobertura
   farmacológica eram indistinguíveis, no v1, de um moderado apoiado nas duas camadas.

## 5. Limites

1. **Isto não melhora o acerto, e não foi validado por acerto.** É uma mudança de relatório: torna
   auditável o que o nível único escondia. O painel tem 15 positivas e 3 negativas avaliadas, com a
   classe atribuída por molécula e não por receptor, então não serve para comparar v1 e v2 por
   sensibilidade. Para isso faltam rótulos externos por receptor, que o eixo de triagem não tem
   (o rótulo do Tox21 é o próprio sinal de treino).
2. **Os limiares de probabilidade continuam os de cada modelo** (0,11 a 0,40 no Tox21 contra 0,5 no
   ChEMBL). A calibração é o passo 2 da proposta e não está implementada; até lá, a probabilidade de
   um eixo não é comparável com a do outro.
3. **A direção vem só do Tox21.** O eixo farmacológico fica como `nao determinada`, porque a anotação
   de modo do ChEMBL é inferida por expressão regular e não é confiável nos ensaios de ligação.
4. **O v2 não produz um nível único e não ordena moléculas.** Se você precisar de uma fila de
   priorização, essa regra de ordenação é uma decisão separada, que eu não tomei aqui.
5. **O painel é de químicos ambientais**, quase todos fora do domínio ChEMBL. A proporção de
   `nao_avaliado` no eixo farmacológico seria menor com moléculas tipo fármaco. Não medi isso.
6. **A potência não é uma escala única.** pIC50, pKi e pEC50 não são intercambiáveis. Em 15 dos 27
   pares com potência no eixo farmacológico há mais de uma medida disponível, então o `pAct_melhor`
   é o máximo entre escalas diferentes e só faz sentido lido junto de `pAct_melhor_medida`. Use a
   coluna `potencia_por_medida`, que traz a potência de cada medida separadamente.
7. **Nos 504 endpoints implantáveis deste painel, a decisão veio sempre da classificação
   quando havia cobertura** (257 casos) ou de nenhuma
   tarefa (247 fora do domínio). O caminho alternativo
   (endpoint com regressão implantável e sem classificador) existe no código, mas não foi exercitado
   aqui, então não está testado.

## 6. Como usar

```
# 1. previsão com o pipeline original, de dentro da pasta do pipeline (os caminhos padrão importam:
#    predict_edc2.py identifica a camada pelo nome da pasta, então --model-dirs absolutos quebram o rótulo)
cd ~/Downloads/pipeline_edc_qsar
python predict_edc2.py --input painel_ref.csv --name-col nome --out /caminho/saida/painel

# 2. reagregação em dois eixos
python consenso_v2.py /caminho/saida/painel_detalhado.csv /caminho/saida/consenso_v2
```

Saídas: `_endpoints.csv` (uma linha por molécula × receptor × endpoint), `_resumo.csv` (uma linha por
molécula × receptor, com os dois eixos e a cobertura) e `_cobertura.csv` (o teto de cada receptor).
