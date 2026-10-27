#!/usr/bin/env python
"""Preditor de consenso multi-camada para desregulação endócrina.

Camadas independentes, nunca fundidas:
  - ChEMBL IC50   (models_final/)   regressão de pIC50 + classificação
  - ChEMBL Ki/EC50 (models_layers/) regressão de pKi / pEC50 + classificação
  - Tox21 qHTS     (models_tox21/)  classificação por ensaio, química ambiental

Uso:
    python predict_edc2.py --smiles "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1"
    python predict_edc2.py --input moleculas.csv --smiles-col smiles --name-col nome
"""
import argparse, glob, os, sys, json
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import joblib
import curation
import features as F

RECEPTOR = {"AR": "AR", "ESR1": "ERalfa", "ESR2": "ERbeta", "PPARG": "PPARgama",
            "THRA": "THRalfa", "THRB": "THRbeta", "TR": "THR (nao resolve alfa/beta)"}
ENSAIO = {"B": "ligacao", "F": "funcional"}


# Load every .joblib model and tag it with its layer, receptor and task.
# Deployability verdicts from criterio_implantacao.csv, keyed by model id. That table is the
# auditable source (see criterio_implantacao.py); the flag baked into each .joblib came from a
# script that is not in the repository. Falls back to the baked flag when the table is absent.
def _criterio(raiz=None):
    raiz = raiz or os.path.dirname(os.path.abspath(__file__))
    p = os.path.join(raiz, "criterio_implantacao.csv")
    if not os.path.exists(p):
        return {}
    t = pd.read_csv(p)
    return {r.modelo_id: (bool(r.nao_implantavel), str(r.motivo)) for r in t.itertuples()}


# Load every .joblib, tagging each with its layer, receptor and deployability verdict.
def carregar(model_dirs=("models_final", "models_layers", "models_tox21")):
    crit = _criterio()
    camadas = []
    for d in model_dirs:
        # The layer is identified by the FOLDER NAME, and each layer has its own filename
        # grammar. Comparing the path as given would mislabel every absolute path, silently
        # parsing Tox21 models as ChEMBL ones, so compare the basename and refuse unknown
        # folders instead of guessing.
        pasta = os.path.basename(os.path.normpath(d))
        if pasta not in ("models_final", "models_layers", "models_tox21"):
            raise SystemExit(
                f"pasta de modelos nao reconhecida: {d!r} (nome base {pasta!r}).\n"
                "Esperado models_final, models_layers ou models_tox21 - o nome define como o "
                "nome do arquivo e interpretado. Renomeie a pasta ou use um atalho simbolico.")
        arquivos = sorted(glob.glob(os.path.join(d, "*.joblib")))
        if not arquivos:
            print(f"aviso: nenhum .joblib em {d}", file=sys.stderr)
        for p in arquivos:
            m = joblib.load(p)
            base = os.path.basename(p)[:-7]
            rep, motivo = crit.get(base, (bool(m.get("nao_implantavel", False)), ""))
            if pasta == "models_tox21":
                _, rec, modo, cel = base.split("_", 3)
                camadas.append({"camada": "Tox21", "receptor": rec, "detalhe": f"{modo}/{cel}",
                                "tarefa": "classificacao", "modelo": m, "chave": base,
                                "reprovado": rep, "motivo_reprovacao": motivo})
            else:
                partes = base.split("_")
                tarefa = partes[-1]
                if pasta == "models_final":
                    rec, atype, medida = partes[0], partes[1], "IC50"
                else:
                    rec, medida, atype = partes[0], partes[1], partes[2]
                camadas.append({"camada": f"ChEMBL {medida}", "receptor": rec,
                                "detalhe": ENSAIO.get(atype, atype),
                                "tarefa": "regressao" if tarefa == "regressao" else "classificacao",
                                "modelo": m, "chave": base,
                                "reprovado": rep, "motivo_reprovacao": motivo})
    return camadas


# Per-compound difficulty used to scale the normalised conformal interval.
def _dificuldade(est, X):
    sub = getattr(est, "estimators_", None)
    if sub is None:
        return np.ones(len(X))
    try:
        return np.vstack([e.predict(X) for e in np.ravel(sub)]).std(axis=0) + 1e-6
    except Exception:
        return np.ones(len(X))


# Score the molecules with all models and aggregate the per-model votes for each receptor.
def prever(smiles_list, nomes=None, model_dirs=("models_final", "models_layers", "models_tox21")):
    camadas = carregar(model_dirs)
    if not camadas:
        raise SystemExit("nenhum modelo encontrado")
    std = [curation.standardize(s) for s in smiles_list]
    base = pd.DataFrame({"nome": nomes if nomes is not None else smiles_list,
                         "smiles_entrada": smiles_list,
                         "std_smiles": [d["std_smiles"] for d in std],
                         "inchikey": [d["inchikey"] for d in std],
                         "motivo_rejeicao": [d["reject_reason"] for d in std]})
    # The summary groups by (nome, receptor), so two rows sharing a name would be merged
    # into one, mixing the votes of different molecules without any error. Refuse instead.
    # Without --name-col the name IS the SMILES, so a repeated SMILES has the same effect.
    dup = base.loc[base["nome"].duplicated(), "nome"].unique()
    if len(dup):
        amostra = ", ".join(map(str, dup[:5])) + (" ..." if len(dup) > 5 else "")
        raise SystemExit(
            f"{len(dup)} nome(s) repetido(s) na entrada: {amostra}\n"
            "O resumo agrupa por (nome, receptor), entao linhas com o mesmo nome seriam "
            "fundidas numa so e os votos de moleculas diferentes se misturariam.\n"
            "Use --name-col com rotulos unicos, ou remova as linhas duplicadas.")
    val = base["std_smiles"].notna().to_numpy()
    if val.sum() == 0:
        return base, pd.DataFrame()
    X, bits, _, ok = F.featurize(base.loc[val, "std_smiles"].tolist())
    idx = base.index[val][ok]

    linhas = []
    for c in camadas:
        m = c["modelo"]
        sim = F.tanimoto_max_to_train(bits, m["train_bits"])
        dentro = sim >= m["ad_threshold"]
        d = pd.DataFrame({"idx": idx, "camada": c["camada"], "receptor": c["receptor"],
                          "detalhe": c["detalhe"], "tarefa": c["tarefa"],
                          "endpoint": c["chave"], "tanimoto_viz": sim.round(3),
                          "no_dominio": dentro,
                          "nao_implantavel": c["reprovado"],
                          "motivo_reprovacao": c["motivo_reprovacao"]})
        if c["tarefa"] == "regressao":
            yh = m["model"].predict(X)
            d["valor_previsto"] = yh.round(2)
            q = m.get("conformal_q_abs")
            if q is not None:
                if m.get("conformal_modo") == "normalizado":
                    q = m["conformal_q_norm"] * _dificuldade(m["model"], X)
                d["intervalo_inf"] = (yh - q).round(2)
                d["intervalo_sup"] = (yh + q).round(2)
                d["cobertura_scaffold"] = m.get("cobertura_scaffold")
            d["ativo"] = d["valor_previsto"] >= 6.0
        else:
            p = m["model"].predict_proba(X)[:, 1]
            thr = float(m.get("limiar", 0.5))
            d["prob_ativo"] = p.round(3)
            d["limiar"] = thr
            d["ativo"] = p >= thr
        linhas.append(d)
    det = pd.concat(linhas, ignore_index=True)
    det = det.merge(base[["nome", "smiles_entrada", "inchikey"]].reset_index()
                    .rename(columns={"index": "idx"}), on="idx", how="left")

    # --- unidades de evidencia ------------------------------------------------------
    # One model is NOT one piece of evidence. Two kinds of redundancy were measured in this
    # project: regression and classification of one ChEMBL endpoint are the SAME compounds
    # and features differing only by binarisation at pAct >= 6; and the Tox21 assays of one
    # receptor+mode are the same compound library (up to 99.8% shared) run in different cell
    # lines. Counting raw model votes therefore overstates the evidence, so votes are first
    # collapsed into units: one unit per ChEMBL endpoint (receptor x measure x assay type)
    # and one unit per Tox21 receptor+mode.
    det["conta"] = det["no_dominio"] & ~det["nao_implantavel"]
    det["motivo_nao_voto"] = np.where(
        det["nao_implantavel"],
        "modelo reprovado (" + det["motivo_reprovacao"].fillna("sem motivo registrado") + ")",
        np.where(~det["no_dominio"], "molecula fora do dominio de aplicabilidade", ""))
    det["eixo"] = np.where(det["camada"].str.startswith("ChEMBL"), "farmacologico", "triagem")
    det["unidade"] = np.where(
        det["eixo"].eq("farmacologico"), det["endpoint"].str.rsplit("_", n=1).str[0],
        det["receptor"] + " " + det["detalhe"].str.split("/").str[0])

    ev = det[det["conta"]].copy()
    ev["eh_cls"] = ev["tarefa"].eq("classificacao")
    ev["cls_ativo"] = ev["ativo"] & ev["eh_cls"]
    # A unit is positive when: on the pharmacological axis the CLASSIFIER decides (regression
    # supplies potency, not the call) and the regression only decides where no classifier is
    # deployable; on the screening axis any cell line calling it active makes the unit positive.
    u = (ev.groupby(["nome", "receptor", "eixo", "unidade"], as_index=False)
           .agg(n_modelos=("ativo", "size"), algum_ativo=("ativo", "any"),
                tem_cls=("eh_cls", "any"), cls_ativo=("cls_ativo", "any")))
    u["positiva"] = np.where(u["eixo"].eq("farmacologico"),
                             np.where(u["tem_cls"], u["cls_ativo"], u["algum_ativo"]),
                             u["algum_ativo"])

    # Pivot the units onto one column per axis, always returning both even when empty.
    def _por_eixo(coluna, func):
        t = u.pivot_table(index=["nome", "receptor"], columns="eixo", values=coluna,
                          aggfunc=func, fill_value=0)
        for e in ("farmacologico", "triagem"):
            if e not in t.columns:
                t[e] = 0
        return t[["farmacologico", "triagem"]].astype(int)

    pos, aval = _por_eixo("positiva", "sum"), _por_eixo("unidade", "nunique")
    # Every molecule x receptor gets a row, including the pairs where NOTHING could be
    # assessed. Grouping only the voting units would drop those pairs from the output
    # entirely - the same "absence disappears silently" defect the ladder exists to avoid,
    # just moved up to the receptor level.
    pares = pd.MultiIndex.from_product(
        [det["nome"].drop_duplicates().tolist(),
         det.loc[~det["nao_implantavel"], "receptor"].drop_duplicates().sort_values().tolist()],
        names=["nome", "receptor"])
    pos, aval = pos.reindex(pares, fill_value=0), aval.reindex(pares, fill_value=0)
    resumo = pd.DataFrame({
        "unidades_farm_positivas": pos["farmacologico"], "unidades_tri_positivas": pos["triagem"],
        "unidades_farm_avaliadas": aval["farmacologico"], "unidades_tri_avaliadas": aval["triagem"],
    }).reset_index()
    resumo["unidades_positivas"] = resumo.unidades_farm_positivas + resumo.unidades_tri_positivas
    resumo["unidades_avaliadas"] = resumo.unidades_farm_avaliadas + resumo.unidades_tri_avaliadas

    # Ceiling: units a receptor could offer at all, identical for every molecule. A unit left
    # out by the domain filter does NOT become a negative; it simply was not assessed, which
    # is why coverage is reported next to the level instead of being folded into it.
    todas = det[~det["nao_implantavel"]].drop_duplicates(["receptor", "eixo", "unidade"])
    teto = (todas.groupby(["receptor", "eixo"]).unidade.nunique().unstack(fill_value=0))
    for e in ("farmacologico", "triagem"):
        if e not in teto.columns:
            teto[e] = 0
    teto = teto.rename(columns={"farmacologico": "teto_farm", "triagem": "teto_tri"}).reset_index()
    resumo = resumo.merge(teto[["receptor", "teto_farm", "teto_tri"]], on="receptor", how="left")
    resumo["unidades_implantaveis"] = resumo.teto_farm + resumo.teto_tri
    resumo["cobertura"] = (resumo.unidades_avaliadas / resumo.unidades_implantaveis).round(2)

    # Agreement between the two axes, reported beside the level instead of altering it.
    f_pos, t_pos = resumo.unidades_farm_positivas > 0, resumo.unidades_tri_positivas > 0
    f_av, t_av = resumo.unidades_farm_avaliadas > 0, resumo.unidades_tri_avaliadas > 0
    resumo["concordancia"] = np.select(
        [f_pos & t_pos, f_pos & t_av & ~t_pos, t_pos & f_av & ~f_pos, f_pos | t_pos,
         f_av & t_av, f_av | t_av],
        ["convergente", "divergente (farmacologico positivo, triagem negativa)",
         "divergente (triagem positiva, farmacologico negativo)", "eixo unico (outro sem cobertura)",
         "negativo nos dois eixos", "negativo (so um eixo avaliado)"],
        default="sem cobertura")

    # Level ladder. "nao avaliado" is a state the previous rule did not have: it separates
    # "tested and negative" from "nobody could test", which the old "sem alerta" merged.
    A, P = resumo.unidades_avaliadas, resumo.unidades_positivas
    resumo["alerta"] = np.select(
        [(A == 0) | ((P == 0) & (A == 1)), P == 0, f_pos & t_pos, P >= 2],
        ["nao avaliado (cobertura insuficiente)", "sem alerta (negativo nas unidades avaliadas)",
         "alto (convergencia entre os dois eixos)", "moderado (2+ unidades independentes)"],
        default="baixo (evidencia isolada)")

    # Direction of effect: only the screening assays carry the mode, so it is reported from
    # the ACTIVE screening votes and stays undetermined when there is none.
    tox_at = ev[ev["eixo"].eq("triagem") & ev["ativo"]].copy()
    if len(tox_at):
        tox_at["modo"] = tox_at["detalhe"].str.split("/").str[0]
        direcao = (tox_at.groupby(["nome", "receptor"]).modo
                   .apply(lambda s: ", ".join(sorted(set(s)))).rename("direcao").reset_index())
        resumo = resumo.merge(direcao, on=["nome", "receptor"], how="left")
    else:
        resumo["direcao"] = np.nan
    resumo["direcao"] = resumo["direcao"].fillna("nao determinada")
    resumo["fracao_positiva"] = (resumo.unidades_positivas
                                 / resumo.unidades_avaliadas.replace(0, np.nan)).round(2)
    # Two different fractions in one string, so each one carries its own label: the same
    # number appearing twice with different meanings ("1/3 unidades (3/6 cobertas)") reads
    # as if 3 of 6 units had been positive.
    resumo["evidencia"] = ("positivas " + resumo.unidades_positivas.astype(str) + "/"
                           + resumo.unidades_avaliadas.astype(str) + "  ·  cobertas "
                           + resumo.unidades_avaliadas.astype(str) + "/"
                           + resumo.unidades_implantaveis.astype(str))
    cols = ["nome", "receptor", "evidencia", "unidades_positivas", "unidades_avaliadas",
            "unidades_implantaveis", "fracao_positiva", "cobertura", "unidades_farm_positivas",
            "unidades_tri_positivas", "unidades_farm_avaliadas", "unidades_tri_avaliadas",
            "concordancia", "direcao", "alerta"]
    return det, resumo[cols]


# Narrow the detail table for writing. The COMPLETE frame is what prever() returns and what
# consenso_v2 needs (its coverage comes from comparing deployable units against assessed
# ones); the narrowing is a presentation choice applied at write time.
#   votantes    linhas que de fato votaram: no dominio e implantavel  (padrao)
#   implantados os modelos implantados, votando ou nao
#   todos       os 61, inclusive os reprovados
def filtrar_detalhe(det, modo="votantes"):
    if modo not in ("votantes", "implantados", "todos"):
        raise SystemExit(f"modo de detalhe invalido: {modo!r} "
                         "(esperado votantes, implantados ou todos)")
    if modo == "todos":
        return det
    d = det[~det["nao_implantavel"]] if modo == "implantados" else det[det["conta"]]
    # colunas que so existem para explicar quem NAO votou ficam constantes depois do
    # filtro, entao saem junto em vez de virar uma coluna vazia no CSV
    sobra = ["motivo_reprovacao"] + (["motivo_nao_voto"] if modo == "votantes" else [])
    return d.drop(columns=[c for c in sobra if c in d.columns]).reset_index(drop=True)


# Command-line entry point of the multi-layer consensus predictor.
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smiles", nargs="+")
    ap.add_argument("--input")
    ap.add_argument("--smiles-col", default="smiles")
    ap.add_argument("--name-col", default=None)
    ap.add_argument("--model-dirs", nargs="+",
                    default=["models_final", "models_layers", "models_tox21"])
    ap.add_argument("--out", default="predicoes_consenso")
    ap.add_argument("--detalhe", choices=("votantes", "implantados", "todos"),
                    default="votantes",
                    help="quais modelos vao para o _detalhado.csv (padrao: votantes)")
    a = ap.parse_args()
    if a.input:
        df = pd.read_csv(a.input)
        smi = df[a.smiles_col].astype(str).tolist()
        nomes = df[a.name_col].tolist() if a.name_col else None
    elif a.smiles:
        smi, nomes = a.smiles, None
    else:
        ap.error("informe --smiles ou --input")
    det, res = prever(smi, nomes, tuple(a.model_dirs))
    det = filtrar_detalhe(det, a.detalhe)
    det.to_csv(f"{a.out}_detalhado.csv", index=False)
    res.to_csv(f"{a.out}_resumo.csv", index=False)
    print(f"{len(smi)} moléculas -> {a.out}_detalhado.csv / {a.out}_resumo.csv")
    if len(res):
        piv = res.pivot(index="nome", columns="receptor", values="alerta")
        print(piv.fillna("-").to_string())


if __name__ == "__main__":
    main()
