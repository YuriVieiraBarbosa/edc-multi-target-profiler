# Pipeline EDC-QSAR — pacote v2

Este pacote é um **incremento** sobre o seu `pipeline_edc_qsar`, não uma substituição. Ele acrescenta
uma camada de relatório (consenso v2), os scripts de análise e os documentos produzidos nas últimas
sessões. Os 17 scripts na raiz são **cópias byte a byte** dos seus originais, conferidas por SHA-256.

## O que NÃO mudou

Antes de qualquer coisa, para evitar mal-entendido:

- **Nenhum modelo foi retreinado.** Os 61 modelos (18 IC50 + 32 Ki/Kd/EC50 + 11 Tox21) são os mesmos.
- **`predict_edc2.py` não foi alterado.** Ele continua produzindo o `_resumo.csv` com o nível único
  (alto / moderado / baixo / sem alerta), exatamente como antes.
- **Nenhum dataset foi alterado.** A curadoria, os filtros e as features são os mesmos.
- **A regra do `scaffold_split` continua a original.** Ela foi testada (veja
  `documentos/relatorio_teste_justo_scaffold.md`), mas não foi trocada: a decisão é sua.
- **O pré-processamento continua o mesmo.** O `StandardScaler` segue aplicado ao vetor inteiro de
  2.180 colunas nos 10 modelos que escalam, e o descritor `Ipc` não foi transformado.

O consenso v2 é uma forma alternativa de **ler** os mesmos votos. Nada obriga a usá-lo, e o formato
antigo continua disponível em paralelo.

## Instalação

Copie o conteúdo deste pacote para dentro da sua pasta `pipeline_edc_qsar`, mantendo a estrutura. Os
modelos e os datasets já estão lá e não vêm no pacote (são cerca de 110 MB de `.joblib` e 25 MB de
CSV). Depois de copiar, a pasta deve ficar assim:

```
pipeline_edc_qsar/
  *.py                      17 scripts originais (inalterados)
  requirements.txt
  models_final/             18 modelos IC50          (já existentes)
  models_layers/            32 modelos Ki/Kd/EC50    (já existentes)
  models_tox21/             11 modelos Tox21         (já existentes)
  datasets_*/               datasets curados         (já existentes)
  consenso_v2/              NOVO  camada de relatório de dois eixos
  analises/                 NOVO  scripts e resultados das análises
  figuras/                  NOVO  infográficos e geradores
  documentos/               NOVO  relatórios e guias
  exemplo_painel_ref/       NOVO  saída de exemplo nos dois formatos
```

## Uso

**Formato original**, sem mudança nenhuma:

```bash
cd pipeline_edc_qsar
python predict_edc2.py --input moleculas.csv --name-col nome --out resultados/run1
```

**Formato de dois eixos** (roda o preditor e reagrega em seguida):

```bash
cd pipeline_edc_qsar
python consenso_v2/rodar_consenso_v2.py --input moleculas.csv --name-col nome --out resultados/run1
```

Isso gera os arquivos do formato antigo (`run1_detalhado.csv`, `run1_resumo.csv`) **e** os do novo
(`run1_v2_resumo.csv`, `run1_v2_endpoints.csv`, `run1_v2_cobertura.csv`). Se você já tem um
`_detalhado.csv` de uma execução anterior, pode reagregar sem prever de novo:

```bash
python consenso_v2/consenso_v2.py resultados/run1_detalhado.csv resultados/run1_v2
```

> **Atenção:** rode sempre a partir da raiz do pipeline e com os caminhos de modelo padrão. O
> `predict_edc2.py` identifica a camada de cada modelo pelo **nome da pasta**, então passar
> `--model-dirs` com caminho absoluto faz a camada Tox21 ser rotulada como ChEMBL e corrompe os votos.
> O `rodar_consenso_v2.py` já cuida disso.

## O que o consenso v2 muda

| | `_resumo.csv` (original) | `_v2_resumo.csv` |
|---|---|---|
| Unidade de evidência | o modelo: regressão e classificação do mesmo endpoint contam 2 vezes | o endpoint: 1 evidência |
| Papel de cada tarefa | as duas votam igual | a classificação decide a positividade; a regressão dá a potência com intervalo conformal |
| Agregação | um nível único, somando ChEMBL e Tox21 | dois eixos, nunca somados |
| Dentro do eixo | fração de votos ativos | melhor evidência e quantos endpoints a corroboram |
| Modelo fora do domínio | some e dilui a fração | estado `nao_avaliado`, com a cobertura ao lado |
| Direção do efeito | ausente | rótulo vindo do modo dos ensaios Tox21 |
| Par sem cobertura | sem linha no arquivo | linha com `sem cobertura` |

O detalhamento, a tabela de cobertura por receptor e a comparação com o formato antigo estão em
`documentos/consenso_v2_especificacao.md`.

## Conteúdo

- **`consenso_v2/`** — `consenso_v2.py` (a reagregação) e `rodar_consenso_v2.py` (atalho que encadeia
  o preditor e a reagregação).
- **`analises/`** — `ablacao_camadas_painel.py` mede o que cada camada acrescenta ao alerta;
  `teste_justo_scaffold/` traz o teste da regra do `scaffold_split`; `resultados/` tem as tabelas
  geradas por essas análises no painel de referência.
- **`figuras/`** — os dois infográficos (o do pipeline atual e o da arquitetura v2), seus geradores e
  o JSON com os números da segunda figura.
- **`documentos/`** — `consenso_v2_especificacao.md`, `guia_engenharia_edc_qsar.md`,
  `relatorio_edc_qsar.md`, `relatorio_camadas_consenso.md`,
  `relatorio_teste_justo_scaffold.md` e o `README_original.md` do projeto.
- **`documentos/verificacao_auditoria_codigo.md`** — verificação, item a item, da auditoria de código:
  o que foi confirmado no código, o que precisa de ressalva, e dois diagnósticos rodados nos 38
  endpoints (`analises/resultados/diagnostico_scaffold_split.csv` e `diagnostico_aciclicos_split.csv`).
- **`exemplo_painel_ref/`** — a saída de `painel_ref.csv` nos dois formatos, mais a comparação par a
  par. Serve para conferir se a sua instalação reproduz os mesmos números.

## Decisões ainda em aberto

Nenhuma destas foi aplicada, e todas dependem de você:

1. **Dupla contagem.** O `consenso_v2.py` resolve, mas o `predict_edc2.py` continua contando duas
   vezes. Não há integração automática entre os dois.
2. **Regra do `scaffold_split`.** Medida, não trocada. Pôr os grupos grandes no treino melhora pouco
   (RMSE mediano −0,018) e o ganho vem do volume de dados, não da composição.
3. **Calibração das probabilidades** (passo 2 da arquitetura). Até fazê-la, as probabilidades dos dois
   eixos não são comparáveis entre si: os limiares vão de 0,11 a 0,40 no Tox21 contra 0,5 no ChEMBL.
4. **Pesos por qualidade** (passo 3). Só faz sentido com validação externa, que ainda não existe para
   o eixo de triagem.
5. **Divergência dos CSVs do Tox21** em relação ao `LEIAME_modelos_tox21.csv`. Causa não investigada.
6. **Pré-processamento.** Escalar só os descritores, com `log` ou recorte no `Ipc`, exigiria retreinar.
7. **Receptor de insulina (IR).** Citado no escopo do projeto, ainda sem camada.
8. **Moléculas sem scaffold de Murcko** caem todas num pseudo-grupo só. No Tox21 isso é 18% a 23% de
   cada ensaio, e deixa dois ensaios com teste externo degenerado. Veja a seção 4 da verificação.

## Limites que acompanham o pacote

- O consenso v2 **não melhora o acerto e não foi validado por acerto**. Ele torna auditável o que o
  nível único escondia.
- O `painel_ref.csv` (15 positivas e 3 negativas avaliadas, com a classe atribuída por molécula e não
  por receptor) é um teste de fumaça, não uma validação.
- Não há verdade independente para o eixo de triagem: o rótulo do Tox21 é o próprio sinal de treino.
  O BindingDB valida só o eixo farmacológico.
- `pIC50`, `pKi` e `pEC50` não são a mesma escala. No `_v2_resumo.csv`, leia a potência pela coluna
  `potencia_por_medida`, que separa as medidas, e não só pelo `pAct_melhor`.

## Ambiente

Python com RDKit, scikit-learn e XGBoost, conforme `requirements.txt`. As análises foram rodadas em um
ambiente com RDKit 2026.03.6.
