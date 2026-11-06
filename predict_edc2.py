#!/usr/bin/env python
"""Multi-layer consensus predictor for endocrine disruption.

Independent layers, never merged:
  - ChEMBL IC50   (models_final/)   pIC50 regression + classification
  - ChEMBL Ki/EC50 (models_layers/) pKi / pEC50 regression + classification
  - Tox21 qHTS     (models_tox21/)  classification per assay, environmental chemistry

Receptors covered: AR, ERalpha, ERbeta, PPARgamma, THRalpha, THRbeta, and PXR. PXR was added
after the others (pipeline_pxr.py and pipeline_pxr_tox21.py) and has both axes, but only in
agonist mode in Tox21 — there is no PXR antagonism counter-assay in the program.

Usage:
    python predict_edc2.py --smiles "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1"
    python predict_edc2.py --input molecules.csv --smiles-col smiles --name-col name
"""
import argparse, glob, os, sys, json
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import joblib
import curation
import features as F

RECEPTOR = {"AR": "AR", "ESR1": "ERalpha", "ESR2": "ERbeta", "PPARG": "PPARgamma",
            "THRA": "THRalpha", "THRB": "THRbeta", "TR": "THR (does not resolve alpha/beta)",
            "PXR": "PXR (pregnane X receptor)"}
# ChEMBL assay kind suffix. ATTENTION: the assay_type field poorly separates binding from
# functional - there are cell reporter assays labeled 'B' across all targets, and in PXR the
# largest functional group comes labeled 'A'. The label below is that of the field, not a guarantee of the
# what the assay measures; the description of each assay is in the curated dataset.
ASSAY = {"B": "binding", "F": "functional", "A": "ADMET/reporter", "T": "transactivation"}


# Load every .joblib model and tag it with its layer, receptor and task.
# Deployability verdicts from deployment_criterion.csv, keyed by model id. That table is the
# auditable source (see deployment_criterion.py); the flag baked into each .joblib came from a
# script that is not in the repository. Falls back to the baked flag when the table is absent.
def _criterion(root=None):
    root = root or os.path.dirname(os.path.abspath(__file__))
    p = os.path.join(root, "deployment_criterion.csv")
    if not os.path.exists(p):
        return {}
    t = pd.read_csv(p)
    return {r.model_id: (bool(r.not_deployable), str(r.reason)) for r in t.itertuples()}


# Load every .joblib, tagging each with its layer, receptor and deployability verdict.
# A model the criterion rejects is NOT loaded: it never votes, so scoring it would only cost
# time. They live in models_failed/ (documentation only; see archive_failed_models.py). If one
# turns up in a run folder anyway, for example after retraining, it is skipped with a warning.
def load(model_dirs=("models_final", "models_layers", "models_tox21")):
    crit = _criterion()
    layers = []
    ignored = []
    for d in model_dirs:
        # The layer is identified by the FOLDER NAME, and each layer has its own filename
        # grammar. Comparing the path as given would mislabel every absolute path, silently
        # parsing Tox21 models as ChEMBL ones, so compare the basename and refuse unknown
        # folders instead of guessing.
        folder = os.path.basename(os.path.normpath(d))
        if folder not in ("models_final", "models_layers", "models_tox21"):
            raise SystemExit(
                f"model folder not recognized: {d!r} (base name {folder!r}).\n"
                "Expected models_final, models_layers or models_tox21 - the name defines how the "
                "file name is interpreted. Rename the folder or use a symbolic link.")
        files = sorted(glob.glob(os.path.join(d, "*.joblib")))
        if not files:
            print(f"warning: no .joblib in {d}", file=sys.stderr)
        for p in files:
            m = joblib.load(p)
            base = os.path.basename(p)[:-7]
            # the keys of the saved .joblib artifacts keep their original names (limiar, modelo, nao_implantavel, ...) so the existing model files stay readable
            rep, reason = crit.get(base, (bool(m.get("nao_implantavel", False)), ""))
            if rep:
                ignored.append(base)
                continue
            if folder == "models_tox21":
                _, rec, mode, cel = base.split("_", 3)
                layers.append({"layer": "Tox21", "receptor": rec, "detail": f"{mode}/{cel}",
                                "task": "classification", "model_name": m, "key": base,
                                "failed": rep, "failure_reason": reason})
            else:
                parts = base.split("_")
                task = parts[-1]
                if folder == "models_final":
                    rec, atype, measure = parts[0], parts[1], "IC50"
                else:
                    rec, measure, atype = parts[0], parts[1], parts[2]
                layers.append({"layer": f"ChEMBL {measure}", "receptor": rec,
                                "detail": ASSAY.get(atype, atype),
                                "task": "regression" if task == "regression" else "classification",
                                "model_name": m, "key": base,
                                "failed": rep, "failure_reason": reason})
    if ignored:
        print(f"warning: {len(ignored)} failed model(s) were in the running folders and "
              f"were ignored ({', '.join(ignored[:3])}{', ...' if len(ignored) > 3 else ''}). "
              "Run python archive_failed_models.py --apply.", file=sys.stderr)
    return layers


# Per-compound difficulty used to scale the normalised conformal interval.
def _difficulty(est, X):
    sub = getattr(est, "estimators_", None)
    if sub is None:
        return np.ones(len(X))
    try:
        return np.vstack([e.predict(X) for e in np.ravel(sub)]).std(axis=0) + 1e-6
    except Exception:
        return np.ones(len(X))


# Score the molecules with all models and aggregate the per-model votes for each receptor.
def predict(smiles_list, names=None, model_dirs=("models_final", "models_layers", "models_tox21")):
    layers = load(model_dirs)
    if not layers:
        raise SystemExit("no model found")
    std = [curation.standardize(s) for s in smiles_list]
    base = pd.DataFrame({"name": names if names is not None else smiles_list,
                         "input_smiles": smiles_list,
                         "std_smiles": [d["std_smiles"] for d in std],
                         "inchikey": [d["inchikey"] for d in std],
                         "rejection_reason": [d["reject_reason"] for d in std]})
    # The summary groups by (name, receptor), so two rows sharing a name would be merged
    # into one, mixing the votes of different molecules without any error. Refuse instead.
    # Without --name-col the name IS the SMILES, so a repeated SMILES has the same effect.
    dup = base.loc[base["name"].duplicated(), "name"].unique()
    if len(dup):
        sample = ", ".join(map(str, dup[:5])) + (" ..." if len(dup) > 5 else "")
        raise SystemExit(
            f"{len(dup)} repeated name(s) in input: {sample}\n"
            "The summary groups by (name, receptor), so rows with the same name would be "
            "merged into one and the votes of different molecules would get mixed together.\n"
            "Use --name-col with unique labels, or remove the duplicate rows.")
    val = base["std_smiles"].notna().to_numpy()
    if val.sum() == 0:
        return base, pd.DataFrame()
    X, bits, _, ok = F.featurize(base.loc[val, "std_smiles"].tolist())
    idx = base.index[val][ok]

    rows = []
    for c in layers:
        m = c["model_name"]
        sim = F.tanimoto_max_to_train(bits, m["train_bits"])
        inside = sim >= m["ad_threshold"]
        d = pd.DataFrame({"idx": idx, "layer": c["layer"], "receptor": c["receptor"],
                          "detail": c["detail"], "task": c["task"],
                          "endpoint": c["key"], "tanimoto_nn": sim.round(3),
                          "in_domain": inside,
                          "not_deployable": c["failed"],
                          "failure_reason": c["failure_reason"]})
        if c["task"] == "regression":
            yh = m["model"].predict(X)
            d["predicted_value"] = yh.round(2)
            # the keys of the saved .joblib artifacts keep their original names (limiar, modelo, nao_implantavel, ...) so the existing model files stay readable
            q = m.get("conformal_q_abs")
            if q is not None:
                if m.get("conformal_modo") == "normalizado":
                    q = m["conformal_q_norm"] * _difficulty(m["model"], X)
                d["interval_lower"] = (yh - q).round(2)
                d["interval_upper"] = (yh + q).round(2)
                d["coverage_scaffold"] = m.get("cobertura_scaffold")
            d["active"] = d["predicted_value"] >= 6.0
        else:
            p = m["model"].predict_proba(X)[:, 1]
            thr = float(m.get("limiar", 0.5))
            d["prob_active"] = p.round(3)
            d["threshold"] = thr
            d["active"] = p >= thr
        rows.append(d)
    det = pd.concat(rows, ignore_index=True)
    det = det.merge(base[["name", "input_smiles", "inchikey"]].reset_index()
                    .rename(columns={"index": "idx"}), on="idx", how="left")

    # --- evidence units ------------------------------------------------------
    # One model is NOT one piece of evidence. Two kinds of redundancy were measured in this
    # project: regression and classification of one ChEMBL endpoint are the SAME compounds
    # and features differing only by binarisation at pAct >= 6; and the Tox21 assays of one
    # receptor+mode are the same compound library (up to 99.8% shared) run in different cell
    # lines. Counting raw model votes therefore overstates the evidence, so votes are first
    # collapsed into units: one unit per ChEMBL endpoint (receptor x measure x assay type)
    # and one unit per Tox21 receptor+mode.
    det["counted"] = det["in_domain"] & ~det["not_deployable"]
    det["reason_not_voted"] = np.where(
        det["not_deployable"],
        "model failed (" + det["failure_reason"].fillna("no reason recorded") + ")",
        np.where(~det["in_domain"], "molecule outside the applicability domain", ""))
    det["axis"] = np.where(det["layer"].str.startswith("ChEMBL"), "pharmacological", "screening")
    det["unit"] = np.where(
        det["axis"].eq("pharmacological"), det["endpoint"].str.rsplit("_", n=1).str[0],
        det["receptor"] + " " + det["detail"].str.split("/").str[0])

    ev = det[det["counted"]].copy()
    ev["is_cls"] = ev["task"].eq("classification")
    ev["cls_active"] = ev["active"] & ev["is_cls"]
    # A unit is positive when: on the pharmacological axis the CLASSIFIER decides (regression
    # supplies potency, not the call) and the regression only decides where no classifier is
    # deployable; on the screening axis any cell line calling it active makes the unit positive.
    u = (ev.groupby(["name", "receptor", "axis", "unit"], as_index=False)
           .agg(n_models=("active", "size"), any_active=("active", "any"),
                has_cls=("is_cls", "any"), cls_active=("cls_active", "any")))
    u["positive"] = np.where(u["axis"].eq("pharmacological"),
                             np.where(u["has_cls"], u["cls_active"], u["any_active"]),
                             u["any_active"])

    # Pivot the units onto one column per axis, always returning both even when empty.
    def _by_axis(column, func):
        t = u.pivot_table(index=["name", "receptor"], columns="axis", values=column,
                          aggfunc=func, fill_value=0)
        for e in ("pharmacological", "screening"):
            if e not in t.columns:
                t[e] = 0
        return t[["pharmacological", "screening"]].astype(int)

    pos, n_assessed = _by_axis("positive", "sum"), _by_axis("unit", "nunique")
    # Every molecule x receptor gets a row, including the pairs where NOTHING could be
    # assessed. Grouping only the voting units would drop those pairs from the output
    # entirely - the same "absence disappears silently" defect the ladder exists to avoid,
    # just moved up to the receptor level.
    pairs = pd.MultiIndex.from_product(
        [det["name"].drop_duplicates().tolist(),
         det.loc[~det["not_deployable"], "receptor"].drop_duplicates().sort_values().tolist()],
        names=["name", "receptor"])
    pos, n_assessed = pos.reindex(pairs, fill_value=0), n_assessed.reindex(pairs, fill_value=0)
    summary = pd.DataFrame({
        "pharm_units_positive": pos["pharmacological"], "screen_units_positive": pos["screening"],
        "pharm_units_assessed": n_assessed["pharmacological"], "screen_units_assessed": n_assessed["screening"],
    }).reset_index()
    summary["units_positive"] = summary.pharm_units_positive + summary.screen_units_positive
    summary["units_assessed"] = summary.pharm_units_assessed + summary.screen_units_assessed

    # Ceiling: units a receptor could offer at all, identical for every molecule. A unit left
    # out by the domain filter does NOT become a negative; it simply was not assessed, which
    # is why coverage is reported next to the level instead of being folded into it.
    all_items = det[~det["not_deployable"]].drop_duplicates(["receptor", "axis", "unit"])
    cap = (all_items.groupby(["receptor", "axis"]).unit.nunique().unstack(fill_value=0))
    for e in ("pharmacological", "screening"):
        if e not in cap.columns:
            cap[e] = 0
    cap = cap.rename(columns={"pharmacological": "cap_pharm", "screening": "cap_screen"}).reset_index()
    summary = summary.merge(cap[["receptor", "cap_pharm", "cap_screen"]], on="receptor", how="left")
    summary["units_deployable"] = summary.cap_pharm + summary.cap_screen
    summary["coverage"] = (summary.units_assessed / summary.units_deployable).round(2)

    # Agreement between the two axes, reported beside the level instead of altering it.
    f_pos, t_pos = summary.pharm_units_positive > 0, summary.screen_units_positive > 0
    f_av, t_av = summary.pharm_units_assessed > 0, summary.screen_units_assessed > 0
    summary["agreement"] = np.select(
        [f_pos & t_pos, f_pos & t_av & ~t_pos, t_pos & f_av & ~f_pos, f_pos | t_pos,
         f_av & t_av, f_av | t_av],
        ["convergent", "divergent (pharmacological positive, screening negative)",
         "divergent (screening positive, pharmacological negative)", "single axis (other without coverage)",
         "negative on both axes", "negative (only one axis assessed)"],
        default="no coverage")

    # Level ladder. "not assessed" is a state the previous rule did not have: it separates
    # "tested and negative" from "nobody could test", which the old "no alert" merged.
    A, P = summary.units_assessed, summary.units_positive
    summary["alert"] = np.select(
        [(A == 0) | ((P == 0) & (A == 1)), P == 0, f_pos & t_pos, P >= 2],
        ["not assessed (insufficient coverage)", "no alert (negative in the units assessed)",
         "high (convergence between the two axes)", "moderate (2+ independent units)"],
        default="low (isolated evidence)")

    # Direction of effect: only the screening assays carry the mode, so it is reported from
    # the ACTIVE screening votes and stays undetermined when there is none.
    tox_at = ev[ev["axis"].eq("screening") & ev["active"]].copy()
    if len(tox_at):
        tox_at["mode"] = tox_at["detail"].str.split("/").str[0]
        direction = (tox_at.groupby(["name", "receptor"])["mode"]
                   .apply(lambda s: ", ".join(sorted(set(s)))).rename("direction").reset_index())
        summary = summary.merge(direction, on=["name", "receptor"], how="left")
    else:
        summary["direction"] = np.nan
    summary["direction"] = summary["direction"].fillna("undetermined")
    summary["positive_fraction"] = (summary.units_positive
                                 / summary.units_assessed.replace(0, np.nan)).round(2)
    # Only what was actually assessed. Coverage (how many units existed and how many
    # could weigh in) stays in the numeric columns units_deployable and coverage, which
    # remain in the CSV but are dropped from the printed table.
    summary["evidence"] = ("positive " + summary.units_positive.astype(str) + "/"
                           + summary.units_assessed.astype(str))
    cols = ["name", "receptor", "evidence", "units_positive", "units_assessed",
            "units_deployable", "positive_fraction", "coverage", "pharm_units_positive",
            "screen_units_positive", "pharm_units_assessed", "screen_units_assessed",
            "agreement", "direction", "alert"]
    return det, summary[cols]


# Narrow the detail table for writing. The COMPLETE frame is what predict() returns and what
# consensus_v2 needs (its coverage comes from comparing deployable units against assessed
# ones); the narrowing is a presentation choice applied at write time.
#   voters    rows that actually voted: in-domain and deployable  (default)
#   deployed all loaded models, voting or not
#   all       same as deployed: the failed ones are no longer loaded
#               (they stay in models_failed/), so there is nothing to add. The name
#               remains accepted so as not to break those who already use it.
def filter_detail(det, mode="voters"):
    if mode not in ("voters", "deployed", "all"):
        raise SystemExit(f"invalid detail mode: {mode!r} "
                         "(expected voters, deployed or all)")
    if mode == "all":
        return det
    d = det[~det["not_deployable"]] if mode == "deployed" else det[det["counted"]]
    # columns that only exist to explain who did NOT vote stay constant after the
    # filter, so they are dropped along with it instead of becoming an empty column in the CSV
    leftover_cols = ["failure_reason"] + (["reason_not_voted"] if mode == "voters" else [])
    return d.drop(columns=[c for c in leftover_cols if c in d.columns]).reset_index(drop=True)


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
    ap.add_argument("--out", default="predictions_consensus")
    ap.add_argument("--detail", choices=("voters", "deployed", "all"),
                    default="voters",
                    help="which models go into the _detail.csv: voters (default), or "
                         "deployed/all to also include those left outside the domain")
    a = ap.parse_args()
    if a.input:
        df = pd.read_csv(a.input)
        smi = df[a.smiles_col].astype(str).tolist()
        names = df[a.name_col].tolist() if a.name_col else None
    elif a.smiles:
        smi, names = a.smiles, None
    else:
        ap.error("provide --smiles or --input")
    det, res = predict(smi, names, tuple(a.model_dirs))
    det = filter_detail(det, a.detail)
    det.to_csv(f"{a.out}_detail.csv", index=False)
    res.to_csv(f"{a.out}_summary.csv", index=False)
    print(f"{len(smi)} molecules -> {a.out}_detail.csv / {a.out}_summary.csv")
    if len(res):
        piv = res.pivot(index="name", columns="receptor", values="alert")
        print(piv.fillna("-").to_string())


if __name__ == "__main__":
    main()
