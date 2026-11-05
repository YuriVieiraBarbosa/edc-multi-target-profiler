#!/usr/bin/env python
"""Ponto de entrada do pipeline EDC-QSAR.

Tres comandos:
    testar    prediz uma molecula ou uma lista de moleculas
    explicar  como usar o programa e como ler o resultado
    ajuda     lista os comandos acima

Exemplos:
    python edc.py testar --smiles "CC(=O)Oc1ccccc1C(=O)O"
    python edc.py testar --arquivo moleculas.csv --coluna-nome nome --eixos
    python edc.py explicar
    python edc.py ajuda
"""
import argparse
import datetime
import glob
import os
import sys

import pandas as pd

AQUI = os.path.dirname(os.path.abspath(__file__))
PASTAS_MODELO = ("models_final", "models_layers", "models_tox21")
ROTULO_CAMADA = {"models_final": "ChEMBL IC50", "models_layers": "ChEMBL Ki/Kd/EC50",
                 "models_tox21": "Tox21 (triagem)"}


# Count the .joblib files of each model directory, without loading them.
def _inventario():
    """Por pasta: (treinados, implantados). Treinado nao e implantado - os reprovados
    por criterio_implantacao.csv ficam em models_reprovados/, so como documentacao, e
    nunca sao carregados. 'Treinados' soma as duas pastas; 'implantados' so a de execucao."""
    import pandas as pd
    p = os.path.join(AQUI, "criterio_implantacao.csv")
    rep = set()
    if os.path.exists(p):
        t = pd.read_csv(p)
        rep = set(t.loc[t.nao_implantavel, "modelo_id"])
    inv = {}
    for d in PASTAS_MODELO:
        ids = [os.path.basename(f)[:-7] for f in glob.glob(os.path.join(AQUI, d, "*.joblib"))]
        arq = glob.glob(os.path.join(AQUI, "models_reprovados", d, "*.joblib"))
        inv[d] = (len(ids) + len(arq), sum(1 for i in ids if i not in rep))
    return inv


# Receptors that were trained but have no deployable model, read from the archive's documentation.
def _receptores_sem_modelo(det):
    import pandas as pd
    p = os.path.join(AQUI, "models_reprovados", "LEIAME_reprovados.csv")
    if not os.path.exists(p):
        return []
    # o receptor e o primeiro trecho do id do modelo, a mesma regra de predict_edc2.carregar
    documentados = {m.split("_")[0] for m in pd.read_csv(p).modelo_id}
    return sorted(documentados - set(det.receptor))


# Resolve the output prefix; with no --saida, create resultados/<timestamp>_<name>/.
def _prefixo_saida(saida, nome_padrao="predicoes"):
    if saida:
        destino = os.path.abspath(saida)
        os.makedirs(os.path.dirname(destino) or ".", exist_ok=True)
        return destino
    agora = datetime.datetime.now()
    pasta = os.path.join(os.getcwd(), "resultados", f"{agora:%Y%m%d_%H%M%S}_{nome_padrao}")
    os.makedirs(pasta, exist_ok=True)
    return os.path.join(pasta, nome_padrao)


# Read the molecules from --smiles or from a CSV file, returning SMILES and names.
def _ler_entrada(a, ap):
    if a.smiles:
        # sem --arquivo nao ha coluna de nome; rotular mol_1..mol_n evita uma tabela
        # com o SMILES inteiro em cada linha, e a legenda abaixo preserva o vinculo
        nomes = [f"mol_{i + 1}" for i in range(len(a.smiles))]
        print("Moleculas recebidas:")
        for n, s in zip(nomes, a.smiles):
            print(f"  {n} = {s}")
        return list(a.smiles), nomes
    if a.arquivo:
        if not os.path.exists(a.arquivo):
            ap.error(f"arquivo nao encontrado: {a.arquivo}")
        df = pd.read_csv(a.arquivo)
        if a.coluna_smiles not in df.columns:
            ap.error(f"a coluna '{a.coluna_smiles}' nao existe em {a.arquivo}. "
                     f"Colunas disponiveis: {', '.join(df.columns)}")
        nomes = None
        if a.coluna_nome:
            if a.coluna_nome not in df.columns:
                ap.error(f"a coluna '{a.coluna_nome}' nao existe em {a.arquivo}")
            nomes = df[a.coluna_nome].astype(str).tolist()
        return df[a.coluna_smiles].astype(str).tolist(), nomes
    ap.error("informe --smiles ou --arquivo")


# ANSI colours, by alert level and by two-axis reading.
RESET = "\033[0m"
COR_ALERTA = {"alto": "\033[1;97;41m",      # branco sobre vermelho
              "moderado": "\033[1;33m",      # amarelo
              "baixo": "\033[36m",           # ciano
              "sem alerta": "\033[90m",      # cinza
              "nao avaliado": "\033[35m"}    # magenta: ausencia de evidencia, nao evidencia de ausencia
COR_LEITURA = (("convergente positivo", "\033[1;97;41m"),
               ("divergente", "\033[1;35m"),
               ("positivo", "\033[1;33m"),
               ("sem cobertura", "\033[90m"),
               ("negativo", "\033[90m"))


# Decide whether to emit ANSI colour: never when piped to a file or when the user opted out.
def _cor_ativa(sem_cor):
    if sem_cor or os.environ.get("NO_COLOR"):
        return False
    return bool(getattr(sys.stdout, "isatty", lambda: False)())


# Pick the colour of one cell: exact alert level, or first matching substring of a reading.
def _cor_de(valor, coluna):
    texto = str(valor)
    if coluna == "alerta":
        return COR_ALERTA.get(texto.split(" (")[0], "")
    for chave, cor in COR_LEITURA:
        if chave in texto:
            return cor
    return ""


# Render a table, colouring one column. Widths come from the PLAIN text, and the escape
# codes are added after padding, so the alignment does not count invisible characters.
def _mostrar(df, coluna_cor=None, usar_cor=False, limite=25, arquivo=""):
    vis = df.head(limite)
    cols = list(vis.columns)
    # float com 2 casas; o resto como texto. str(1.0) sairia "1.0" e quebraria a coluna.
    def _fmt(v):
        if isinstance(v, float):
            # v != v e NaN: aparece quando nada foi avaliado (0 de 0), e "nan" na tela
            # parece defeito em vez de ausencia
            return f"{v:.2f}" if v == v else "-"
        return str(v)
    txt = {c: [_fmt(v) for v in vis[c]] for c in cols}
    larg = {c: max([len(c)] + [len(v) for v in txt[c]]) for c in cols}
    print("  ".join(c.rjust(larg[c]) for c in cols))
    for i in range(len(vis)):
        celulas = []
        for c in cols:
            bruto = txt[c][i]
            cor = _cor_de(bruto, c) if (usar_cor and c == coluna_cor) else ""
            # o preenchimento fica FORA do codigo de cor, senao o fundo colorido do
            # nivel "alto" se estende pelos espacos e vira uma barra.
            celulas.append(" " * (larg[c] - len(bruto)) + f"{cor}{bruto}{RESET}" if cor
                           else bruto.rjust(larg[c]))
        print("  ".join(celulas))
    if len(df) > limite:
        print(f"... e mais {len(df) - limite} linha(s). A tabela completa esta em {arquivo}")


# One-line colour legend, printed only when colour is on.
def _legenda_cores(coluna):
    if coluna == "alerta":
        itens = [(n, COR_ALERTA[n]) for n in ("alto", "moderado", "baixo", "sem alerta", "nao avaliado")]
    else:
        itens = [("convergente positivo", COR_LEITURA[0][1]), ("divergente", COR_LEITURA[1][1]),
                 ("so um eixo positivo", COR_LEITURA[2][1]), ("negativo / sem cobertura", COR_LEITURA[3][1])]
    print("legenda: " + "   ".join(f"{cor}{nome}{RESET}" for nome, cor in itens))


# COMMAND 1: predict one molecule or a list of molecules and print the consensus table.
def cmd_testar(a, ap):
    cor = _cor_ativa(a.sem_cor)
    smi, nomes = _ler_entrada(a, ap)
    # predict_edc2.py identifica a camada pelo NOME da pasta de modelos, com caminho relativo;
    # por isso a predicao roda a partir da raiz do pipeline.
    anterior = os.getcwd()
    sys.path.insert(0, AQUI)
    os.chdir(AQUI)
    try:
        import predict_edc2
        det, res = predict_edc2.prever(smi, nomes)
    finally:
        os.chdir(anterior)

    prefixo = _prefixo_saida(a.saida)
    det_completo = det
    det = predict_edc2.filtrar_detalhe(det, a.detalhe)
    det.to_csv(f"{prefixo}_detalhado.csv", index=False)
    res.to_csv(f"{prefixo}_resumo.csv", index=False)

    # Moleculas que nao chegaram a ser pontuadas: a padronizacao estrutural as rejeitou
    # (SMILES invalido, elemento inorganico, sal sem fragmento organico, massa fora de
    # 100-1000 Da). predict_edc2.prever so devolve o motivo quando TODAS falham, entao aqui
    # o motivo e recalculado apenas para as que faltam -- barato, e torna a mensagem util.
    nomes_in = list(nomes) if nomes is not None else list(smi)
    pontuadas = set(det.nome.unique()) if "endpoint" in det.columns else set()
    faltando = [(n, s) for n, s in zip(nomes_in, smi) if n not in pontuadas]
    rej = pd.DataFrame(columns=["nome", "smiles_entrada", "motivo_rejeicao"])
    if faltando:
        sys.path.insert(0, AQUI)
        import curation
        rej = pd.DataFrame([{"nome": n, "smiles_entrada": s,
                             "motivo_rejeicao": curation.standardize(s)["reject_reason"]}
                            for n, s in faltando])
        rej.to_csv(f"{prefixo}_rejeitadas.csv", index=False)

    print(f"\n{len(smi)} molecula(s) enviada(s); {len(pontuadas)} pontuada(s) pelos modelos.")
    if len(rej):
        print(f"{len(rej)} rejeitada(s) na padronizacao estrutural:")
        for r in rej.itertuples():
            print(f"   {r.nome}: {r.motivo_rejeicao}")
        print(f"   (detalhe em {prefixo}_rejeitadas.csv)")
    if len(pontuadas) == 0:
        print("Nenhuma molecula chegou aos modelos. Corrija a entrada e rode de novo.")
    elif len(res) == 0:
        print("Nenhum modelo pode votar: todas as moleculas ficaram fora do dominio "
              "de aplicabilidade de todos os modelos implantaveis.")
    else:
        colunas = ["nome", "receptor", "evidencia", "fracao_positiva",
                   "concordancia", "direcao", "alerta"]
        print()
        votaram = int(res.unidades_avaliadas.sum())
        teto = int(res.unidades_implantaveis.sum())
        impl = int((~det_completo["nao_implantavel"]).sum())
        vot = int(det_completo["conta"].sum())
        print(f"Votaram {vot} de {impl} modelos implantados, que formam {votaram} de {teto} "
              f"unidades de evidencia ({votaram / teto:.0%} de cobertura).")
        if a.detalhe == "votantes":
            print(f"O _detalhado.csv traz so esses {vot}. Use --detalhe implantados "
                  "para ver os que ficaram fora do dominio e por que.")
        _mostrar(res[colunas], coluna_cor="alerta", usar_cor=cor,
                 arquivo=f"{prefixo}_resumo.csv")
        if cor:
            _legenda_cores("alerta")
        sem_cobertura = res[res.unidades_avaliadas == 0]
        if len(sem_cobertura):
            print(f"\nAtencao: {len(sem_cobertura)} par(es) molecula-receptor sem nenhum voto "
                  "(molecula fora do dominio de todos os modelos daquele receptor).")

    if a.eixos:
        caminho_v2 = os.path.join(AQUI, "consenso_v2", "consenso_v2.py")
        if not os.path.exists(caminho_v2):
            print("\n[--eixos ignorado] consenso_v2/consenso_v2.py nao encontrado.")
        else:
            sys.path.insert(0, os.path.dirname(caminho_v2))
            import consenso_v2
            # a tabela COMPLETA, porque a cobertura do v2 compara implantaveis com avaliados
            e = consenso_v2.por_endpoint(consenso_v2.carregar(det_completo), unidade=a.unidade)
            r = consenso_v2.resumir(e, receptores_sem_modelo=_receptores_sem_modelo(det_completo))
            r.to_csv(f"{prefixo}_eixos_resumo.csv", index=False)
            e.to_csv(f"{prefixo}_eixos_unidades.csv", index=False)
            print(f"\nLeitura em dois eixos (unidade de voto: {a.unidade}):")
            _mostrar(r[["nome", "receptor", "leitura", "direcao"]], coluna_cor="leitura",
                     usar_cor=cor, arquivo=f"{prefixo}_eixos_resumo.csv")
            if cor:
                _legenda_cores("leitura")

    print(f"\nArquivos gravados com o prefixo: {prefixo}")
    return 0


TEXTO_EXPLICACAO = """
==============================================================================
 EDC-QSAR - como usar e como ler o resultado
==============================================================================

O QUE O PROGRAMA FAZ
  Recebe uma molecula (SMILES) e devolve, para cada receptor, uma indicacao de
  possivel desregulacao endocrina. Ele NAO mede nada: e uma triagem
  computacional para priorizar compostos, nunca uma prova de atividade.

  Receptores cobertos: AR, ERalfa (ESR1), ERbeta (ESR2), PPARgama, THRalfa,
  THRbeta e PXR (receptor de pregnano X, NR1I2).

  O PXR nao e um receptor hormonal como os outros: e o xenossensor que regula a
  CYP3A4. Ele entra porque medeia parte do efeito dos desreguladores por inducao
  enzimatica e compartilha ligantes com ER, AR e PPARgama. No PXR a triagem so
  tem modo agonista, entao a direcao dele nunca sai como antagonista, e o sufixo
  de tipo de ensaio (A, B, F) nao indica se o ensaio e de ligacao ou funcional.

COMO USAR
  1) Uma molecula:
       python edc.py testar --smiles "CC(=O)Oc1ccccc1C(=O)O"
  2) Varias moleculas de uma vez:
       python edc.py testar --smiles "SMILES_1" "SMILES_2" "SMILES_3"
  3) Um arquivo CSV (recomendado para listas):
       python edc.py testar --arquivo moleculas.csv --coluna-nome nome
     O CSV precisa de uma coluna de SMILES (padrao: 'smiles'). Use
     --coluna-smiles se ela tiver outro nome, e --coluna-nome para o rotulo.
  4) Com a leitura em dois eixos, que separa evidencia farmacologica de
     triagem ambiental:
       python edc.py testar --arquivo moleculas.csv --eixos

COMO O RESULTADO E FORMADO
  Cada modelo implantavel da um voto binario (ativo / inativo):
    - regressao     vota ativo se o valor previsto for >= 6,0 (potencia <= 1 uM)
    - classificacao vota ativo se a probabilidade passar o limiar do modelo
  So vota o modelo que for implantavel E cuja molecula esteja dentro do seu
  dominio de aplicabilidade. Um modelo fora do dominio NAO vota, e tambem NAO
  conta como voto inativo - por isso a cobertura aparece ao lado do nivel.

  O _detalhado.csv traz por padrao SO os modelos que votaram. Para ver tambem os
  que ficaram de fora do dominio e o motivo (coluna motivo_nao_voto), use
  --detalhe implantados. Os modelos reprovados nao sao carregados: ficam em
  models_reprovados/ so como documentacao, com o motivo de cada um.

UNIDADES DE EVIDENCIA - por que um modelo nao e uma evidencia
  Contar um voto por modelo exagera a evidencia, porque ha duas redundancias
  medidas neste projeto:
    - regressao e classificacao do mesmo endpoint usam os MESMOS compostos e as
      mesmas features, diferindo so pela binarizacao em pAct >= 6;
    - os ensaios Tox21 de um mesmo receptor e modo sao a MESMA biblioteca de
      compostos (ate 99,8% em comum) rodada em linhagens celulares diferentes.
  Por isso os votos sao primeiro colapsados em UNIDADES:
    - eixo farmacologico: 1 unidade por endpoint ChEMBL (receptor x medida x
      tipo de ensaio). A classificacao decide se a unidade e positiva; a
      regressao entra como potencia prevista, nao como voto separado.
    - eixo de triagem: 1 unidade por receptor e modo (agonista / antagonista),
      colapsando as linhagens. A unidade e positiva se qualquer linhagem acusar.

OS DOIS EIXOS
  O eixo diz DE ONDE vem a evidencia: ChEMBL (farmacologico, afinidade e
  potencia em ensaio dirigido) ou Tox21 (triagem, qHTS em celula). As duas
  coisas nao se substituem, e por isso o nivel mais alto exige as duas.

COLUNAS DO _resumo.csv
  evidencia                ex. "positivas 1/3": 1 unidade positiva entre as 3
                           que puderam ser avaliadas
  unidades_positivas       unidades de evidencia positivas
  unidades_avaliadas       unidades que puderam ser avaliadas
  unidades_implantaveis    teto do receptor: unidades que existem
  fracao_positiva          positivas / avaliadas
  cobertura                avaliadas / implantaveis  (0 a 1)
  unidades_farm_*          o mesmo, so no eixo farmacologico
  unidades_tri_*           o mesmo, so no eixo de triagem
  concordancia             convergente, divergente, eixo unico ou negativo
  direcao                  agonista / antagonista, dos votos ATIVOS do Tox21
  alerta                   alto, moderado, baixo, sem alerta ou nao avaliado

  A tabela impressa mostra so o que foi avaliado. As colunas 'cobertura' e
  'unidades_implantaveis' nao aparecem na tela, mas estao no _resumo.csv - e
  vale consulta-las antes de comparar duas moleculas: uma fracao positiva de
  1,00 com cobertura 0,17 quer dizer que uma unica unidade de sete opinou, nao
  que a evidencia seja unanime. A linha impressa logo acima da tabela da a
  cobertura agregada da execucao.

NIVEIS DE ALERTA (vale a primeira regra que se aplica)
  nao avaliado  nenhuma unidade avaliada, ou uma so unidade e negativa -
                cobertura insuficiente para afirmar qualquer coisa
  sem alerta    todas as unidades avaliadas sao negativas (2 ou mais)
  alto          os DOIS eixos positivos (convergencia)
  moderado      2 ou mais unidades independentes positivas num eixo so
  baixo         exatamente 1 unidade positiva

  'nao avaliado' e um estado que a regra anterior nao tinha: ela chamava de
  "sem alerta" tanto a molecula testada e negativa quanto aquela que modelo
  nenhum conseguiu avaliar. Sao coisas opostas.

COMO INTERPRETAR - o que o alerta NAO diz
  1. A direcao vem SO do Tox21, que e quem registra o modo do ensaio. No eixo
     farmacologico ela fica "nao determinada", e o nivel em si nao separa
     agonismo de antagonismo.
  2. O teto depende do receptor. "Alto" exige convergencia entre os dois
     eixos, entao so e atingivel onde existem os dois: AR, ERalfa, PPARgama e
     PXR - este ultimo o mais magro dos quatro, com 2 unidades farmacologicas e
     1 de triagem, contra 4 e 4 do ERalfa.
     ERbeta nao tem ensaio Tox21 e para em moderado; THRalfa nao tem nenhum
     modelo implantavel e nao aparece. A coluna unidades_implantaveis mostra
     esse teto, e por isso ela deve ser lida junto do nivel.
  3. "Ativo" nao significa o mesmo em cada camada. No ChEMBL e potencia de
     ate 1 uM; no Tox21 e a chamada ativo/inativo do ensaio, com limiar de
     probabilidade mais permissivo.
  4. Os votos nao sao independentes. Ki e EC50 do mesmo receptor compartilham
     compostos, e os ensaios Tox21 de um receptor sao os mesmos modos em
     linhagens celulares diferentes.
  5. Perto do corte de 6,0 o voto da regressao e ruido. Prefira a
     probabilidade da classificacao como indicador de atividade, e leia a
     regressao como ordem de grandeza da potencia, junto do intervalo
     conformal (colunas intervalo_inf e intervalo_sup do _detalhado.csv).
  6. Cobertura baixa e comum em moleculas tipo farmaco: os modelos do ChEMBL
     sao series estreitas de quimica medicinal, entao muitos farmacos ficam
     fora do dominio deles e o veredito acaba vindo so do Tox21.

ARQUIVOS GERADOS
  <prefixo>_resumo.csv       uma linha por molecula x receptor
  <prefixo>_detalhado.csv    uma linha por modelo, com probabilidade,
                             valor previsto, intervalo conformal, similaridade
                             ao treino e se o voto contou
  com --eixos:
  <prefixo>_eixos_resumo.csv    leitura em dois eixos, sem somar as camadas
  <prefixo>_eixos_unidades.csv  a evidencia por unidade de voto

OUTROS PROGRAMAS DA PASTA
  python rodar_edc.py --input lista.csv --name-col nome --out corrida
      mesma predicao, mais os graficos, arquivada em resultados/
  python consenso_v2/consenso_v2.py <detalhado.csv> <prefixo> [modelo|endpoint]
      so a reagregacao em dois eixos, a partir de um _detalhado.csv existente
==============================================================================
"""


# COMMAND 2: print the usage and interpretation guide, plus a live model inventory.
def cmd_explicar(a, ap):
    print(TEXTO_EXPLICACAO.strip())
    inv = _inventario()
    total = sum(t for t, _ in inv.values())
    impl = sum(i for _, i in inv.values())
    print("\nMODELOS INSTALADOS NESTA PASTA  (implantados / treinados)")
    for d, (t, i) in inv.items():
        print(f"  {ROTULO_CAMADA[d]:22s} {i:3d} / {t:<3d}  ({d}/)")
    print(f"  {'TOTAL':22s} {impl:3d} / {total:<3d}")
    print("  Treinado nao e implantado: os reprovados ficam em models_reprovados/ so como")
    print("  documentacao e nunca sao carregados.")
    print("  Os do Tox21 contam como implantados por AUSENCIA de metrica de validacao")
    print("  gravada, nao por terem passado em um criterio. Motivo por modelo em")
    print("  criterio_implantacao.csv.")
    if total == 0:
        print("\n  ATENCAO: nenhum modelo encontrado. Rode este comando a partir da raiz do"
              "\n  pipeline, ou verifique se as pastas models_* estao no lugar.")
    return 0


# COMMAND 3: list the available commands and the most common options.
def cmd_ajuda(a, ap):
    print("""
EDC-QSAR - comandos disponiveis

  1) testar     prediz uma molecula ou uma lista de moleculas
       python edc.py testar --smiles "C1=CC=CC=C1O"
       python edc.py testar --smiles "SMILES_1" "SMILES_2"
       python edc.py testar --arquivo moleculas.csv --coluna-nome nome
       python edc.py testar --arquivo moleculas.csv --eixos --saida saida/teste1

     opcoes:
       --smiles ...           um ou mais SMILES na linha de comando
       --arquivo CSV          arquivo com uma coluna de SMILES
       --coluna-smiles NOME   nome da coluna de SMILES   (padrao: smiles)
       --coluna-nome NOME     coluna com o rotulo da molecula  (opcional)
       --saida PREFIXO        prefixo dos arquivos de saida
                              (padrao: resultados/<data_hora>_predicoes)
       --eixos                tambem grava a leitura em dois eixos
       --unidade modelo|endpoint   unidade de voto com --eixos (padrao: modelo)
       --sem-cor              desliga as cores (ja desligadas ao redirecionar)
       --detalhe MODO         linhas do _detalhado.csv: votantes (padrao),
                              implantados ou todos

  2) explicar   como usar o programa e como ler o resultado
       python edc.py explicar

  3) ajuda      mostra esta lista
       python edc.py ajuda

Detalhe de cada comando:  python edc.py <comando> --help
""".strip())
    return 0


# Build the argument parser with one subcommand per function.
def construir_parser():
    ap = argparse.ArgumentParser(
        prog="edc.py", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="comando", metavar="comando")

    p1 = sub.add_parser("testar", help="prediz uma molecula ou uma lista de moleculas",
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    p1.add_argument("--smiles", nargs="+", help="um ou mais SMILES")
    p1.add_argument("--arquivo", help="CSV com uma coluna de SMILES")
    p1.add_argument("--coluna-smiles", default="smiles", dest="coluna_smiles",
                    help="nome da coluna de SMILES no CSV (padrao: smiles)")
    p1.add_argument("--coluna-nome", default=None, dest="coluna_nome",
                    help="coluna com o rotulo da molecula (opcional)")
    p1.add_argument("--saida", default=None,
                    help="prefixo dos arquivos de saida (padrao: resultados/<data_hora>_predicoes)")
    p1.add_argument("--eixos", action="store_true",
                    help="tambem grava a leitura em dois eixos (farmacologico x triagem)")
    p1.add_argument("--unidade", choices=("modelo", "endpoint"), default="modelo",
                    help="unidade de voto usada com --eixos (padrao: modelo)")
    p1.add_argument("--detalhe", choices=("votantes", "implantados", "todos"),
                    default="votantes",
                    help="quais modelos vao para o _detalhado.csv (padrao: votantes)")
    p1.add_argument("--sem-cor", action="store_true", dest="sem_cor",
                    help="desliga as cores do terminal (ja desligadas ao redirecionar a saida)")
    p1.set_defaults(func=cmd_testar)

    p2 = sub.add_parser("explicar", help="como usar o programa e como ler o resultado")
    p2.set_defaults(func=cmd_explicar)

    p3 = sub.add_parser("ajuda", help="lista os comandos disponiveis")
    p3.set_defaults(func=cmd_ajuda)
    return ap


# Entry point: with no command, show the command list instead of an error.
def main(argv=None):
    ap = construir_parser()
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        return cmd_ajuda(None, ap)
    a = ap.parse_args(argv)
    if not getattr(a, "func", None):
        return cmd_ajuda(a, ap)
    return a.func(a, ap)


if __name__ == "__main__":
    sys.exit(main())
