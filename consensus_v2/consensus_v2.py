#!/usr/bin/env python
"""Consensus v2: report with two axes, without retraining anything.

Reads the `<out>_detail.csv` produced by predict_edc2.py and rewrites only the aggregation. No model is
altered or retrained, and the individual votes are exactly the same.

Changes relative to the current `_summary.csv`:
  1. The unit of evidence is the ENDPOINT, not the model. The regression and classification of the same ChEMBL
     endpoint no longer count twice: the classification decides positivity (generalizes better outside
     of training) and the regression enters as the order of magnitude of potency, with the conformal interval.
  2. Two axes, NEVER summed: pharmacological (ChEMBL) and screening (Tox21). Each has its own state.
  3. Within the axis, the BEST evidence and how many endpoints corroborate it, instead of the fraction of votes
     (the endpoints are not independent: Ki and EC50 share compounds, and the Tox21 assays of a
     receptor are the same modes in different cell lines).
  4. Three states per axis: positive, negative and `not_assessed`. The coverage is explicit
     (endpoints that exist / that could vote), and whatever is outside the domain never dilutes anything.
  5. The direction of the effect travels as a label (from the mode of the Tox21 assays), outside the score.
  6. Every molecule x receptor combination gets a row, even with no coverage at all. In the current format
     these pairs simply did not appear, and the absence of a row was confused with "no alert".

Usage:
    python consensus_v2.py <out>_detail.csv [output_prefix]

Outputs:
    <prefix>_endpoints.csv  one row per molecule x receptor x endpoint (unit of evidence)
    <prefix>_summary.csv     one row per molecule x receptor, with the two axes and the coverage
    <prefix>_coverage.csv  endpoints available per receptor and axis (the "cap" of each receptor)
"""
import sys

import numpy as np
import pandas as pd

AXIS = {"ChEMBL": "pharmacological", "Tox21": "screening"}
BOOL = ("counted", "active", "in_domain", "not_deployable")


# Read the predictor's _detail.csv and normalise its boolean, optional and derived columns.
def load(path):
    # accepts the CSV path or the DataFrame itself, because the saved detail
    # can be filtered to only the voters and the coverage here needs the whole table
    d = path.copy() if isinstance(path, pd.DataFrame) else pd.read_csv(path)
    # a _detail.csv file written with --detail voters comes already filtered and doesn't include all
    # these columns; in that case the coverage reported here will be that of what remained in the file
    DEFAULT = {"counted": True, "active": False, "in_domain": True, "not_deployable": False}
    for c in BOOL:
        if c not in d.columns:
            d[c] = DEFAULT[c]
        d[c] = d[c].astype(str).str.lower().eq("true")
    # optional columns: predict_edc2.py only creates interval_lower/upper when there is conformal calibration,
    # and prob_active/threshold only when there is a classifier. Without this, a file without conformal breaks below.
    for c in ("predicted_value", "interval_lower", "interval_upper", "prob_active", "threshold", "tanimoto_nn"):
        if c not in d.columns:
            d[c] = np.nan
    if "detail" not in d.columns:
        d["detail"] = ""
    d["axis"] = np.where(d.layer == "Tox21", "screening", "pharmacological")
    # base endpoint: the name without the task suffix (only the ChEMBL layers have both tasks)
    d["endpoint_base"] = d.endpoint.str.replace(r"_(regression|classification)$", "", regex=True)
    d["measure"] = np.where(d.axis == "screening", "AC50", d.layer.str.replace("ChEMBL ", "", regex=False))
    return d


# Build one evidence row per voting unit, for every molecule x receptor x axis.
def by_endpoint(d, unit="model"):
    """One line of evidence per vote unit.

    unit="model"   (default): each model votes alone. Regression and classification of the
                       same endpoint are TWO independent votes. It is the original scheme
                       of predict_edc2.py.
    unit="endpoint": collapses regression + classification of the same endpoint into ONE
                       piece of evidence (the classification decides, the regression becomes potency).
    """
    if unit not in ("model", "endpoint"):
        raise ValueError(f"unit must be 'model' or 'endpoint', received {unit!r}")
    key = ["name", "receptor", "axis", "measure", "endpoint_base"]
    if unit == "model":
        # including the endpoint with task suffix puts each model in its own group
        key = key + ["endpoint"]
    rows = []
    for k, g in d.groupby(key, sort=True):      # sort=True: the order of the rows does not depend on the input
        cls = g[g.task == "classification"]
        reg = g[g.task == "regression"]
        for name_t, t in (("classification", cls), ("regression", reg)):
            if len(t) > 1:                            # never silently choose which row counts
                raise ValueError(
                    f"{len(t)} rows of {name_t} for the same endpoint in {dict(zip(key, k))}. "
                    "Usual cause: two molecules with the same value in the name column (or SMILES "
                    "repeated without --name-col), that predict_edc2.py merges into a single row.")
        is_deployable = bool((~g.not_deployable).any())          # there is some deployable model in this endpoint
        cls_c, reg_c = cls[cls.counted], reg[reg.counted]
        if len(cls_c):                                    # the classification decides positivity
            base, positive = "classification", bool(cls_c.active.iloc[0])
            prob, threshold = float(cls_c.prob_active.iloc[0]), float(cls_c.threshold.iloc[0])
        elif len(reg_c):                                  # no classifier in the domain: falls back to regression
            base, positive = "regression", bool(reg_c.active.iloc[0])
            prob, threshold = np.nan, np.nan
        else:
            base, positive, prob, threshold = "not_assessed", None, np.nan, np.nan
        r = dict(zip(key, k))
        r.update(deployable=is_deployable, assessed=base != "not_assessed", decision_base=base,
                 positive=positive, prob_active=prob, threshold=threshold,
                 tanimoto_nn=float(g.tanimoto_nn.max()),
                 mode=g.detail.iloc[0] if r["axis"] == "screening" else "",
                 pAct_predicted=float(reg_c.predicted_value.iloc[0]) if len(reg_c) else np.nan,
                 pAct_lower=float(reg_c.interval_lower.iloc[0]) if len(reg_c) else np.nan,
                 pAct_upper=float(reg_c.interval_upper.iloc[0]) if len(reg_c) else np.nan)
        rows.append(r)
    e = pd.DataFrame(rows)
    e["simple_mode"] = e["mode"].str.split("/").str[0]
    return e


# State of one axis for a molecule x receptor pair, with vote counts, coverage and best evidence.
def _state(g):
    """State of an axis for a molecule x receptor pair, based on the endpoints of that axis."""
    disp = g[g.deployable]
    evaluated_rows = disp[disp.assessed]
    pos = evaluated_rows[evaluated_rows.positive == True]
    if len(disp) == 0:
        est = "no_model"
    elif len(evaluated_rows) == 0:
        est = "not_assessed"
    else:
        est = "positive" if len(pos) else "negative"
    # best evidence: highest probability among positive endpoints (or among those assessed)
    ref = pos if len(pos) else evaluated_rows
    i = ref.prob_active.idxmax() if len(ref) and ref.prob_active.notna().any() else (ref.index[0] if len(ref) else None)
    best = ref.loc[i] if i is not None else None
    # Potency: pIC50, pKi and pEC50 are NOT on the same scale, so the maximum across measures is not
    # interpretable alone. We store which endpoint and measure it came from, and also the potency of each
    # measure separately, which is the correct reading when the receptor has more than one.
    pot = evaluated_rows[evaluated_rows.pAct_predicted.notna()]
    j = pot.pAct_predicted.idxmax() if len(pot) else None
    by_measure = "; ".join(f"{med} {pot.loc[idx, 'pAct_predicted']:.2f}"
                           for med, idx in sorted(pot.groupby("measure").pAct_predicted.idxmax().items()))
    # "votes received / number of models", in the two possible readings:
    #   positive_fraction           -> denominator = units that voted (previous scheme)
    #   fraction_of_deployable -> denominator = all deployable units of the axis
    # The difference between the two is the coverage of the applicability domain.
    return pd.Series(dict(
        state=est, endpoints_available=len(disp), endpoints_assessed=len(evaluated_rows), endpoints_positive=len(pos),
        positive_fraction=(round(len(pos) / len(evaluated_rows), 2) if len(evaluated_rows) else np.nan),
        fraction_of_deployable=(round(len(pos) / len(disp), 2) if len(disp) else np.nan),
        best_endpoint=(best.endpoint_base if best is not None else ""),
        best_prob=(best.prob_active if best is not None else np.nan),
        best_base=(best.decision_base if best is not None else ""),
        pAct_best=(pot.loc[j, "pAct_predicted"] if j is not None else np.nan),
        pAct_best_lower=(pot.loc[j, "pAct_lower"] if j is not None else np.nan),
        pAct_best_upper=(pot.loc[j, "pAct_upper"] if j is not None else np.nan),
        pAct_best_endpoint=(pot.loc[j, "endpoint_base"] if j is not None else ""),
        pAct_best_measure=(pot.loc[j, "measure"] if j is not None else ""),
        potency_by_measure=by_measure,
        tanimoto_max=(disp.tanimoto_nn.max() if len(disp) else np.nan),
        positive_measures=", ".join(sorted(pos.measure.unique())),
        positive_modes=", ".join(sorted(set(pos.simple_mode.dropna()) - {""})),
    ))


# One row per molecule x receptor, with both axes side by side and the joint reading.
# receptors_without_model: receptors that were trained but have NO deployable model, so no unit of
# theirs reaches `e`. They still get a row (state no_model on both axes, reading "no coverage")
# so that a missing receptor is visible instead of silently absent.
def summarize(e, receptors_without_model=()):
    # indexes once by (name, receptor, axis) instead of refiltering the frame for each pair
    blocks = {k: g for k, g in e.groupby(["name", "receptor", "axis"], sort=False)}
    empty = e.iloc[0:0]
    pairs = e[["name", "receptor"]].drop_duplicates()
    missing_items = [x for x in receptors_without_model if x not in set(e.receptor)]
    if missing_items:
        extra = pd.DataFrame([(mol, x) for mol in pairs.name.unique() for x in missing_items],
                             columns=["name", "receptor"])
        pairs = pd.concat([pairs, extra], ignore_index=True)
    pairs = pairs.sort_values(["name", "receptor"])
    out = []
    for _, (name, rec) in pairs.iterrows():
        r = dict(name=name, receptor=rec)
        for axis, pref in (("pharmacological", "pharm"), ("screening", "tox")):
            s = _state(blocks.get((name, rec, axis), empty))
            r.update({f"{pref}_{k}": v for k, v in s.items()})
        out.append(r)
    r = pd.DataFrame(out)
    # joint reading of the two axes, without summing anything: all 9 state combinations are named
    READING = {
        ("positive", "positive"): "convergent positive",
        ("positive", "negative"): "divergent: pharmacological positive, screening negative",
        ("positive", "no coverage"): "only pharmacological positive (no screening coverage)",
        ("negative", "positive"): "divergent: screening positive, pharmacological negative",
        ("no coverage", "positive"): "only screening positive (no pharmacological coverage)",
        ("negative", "negative"): "convergent negative",
        ("negative", "no coverage"): "negative only on the pharmacological axis (partial coverage)",
        ("no coverage", "negative"): "negative only on the screening axis (partial coverage)",
        ("no coverage", "no coverage"): "no coverage",
    }
    simp = lambda s: np.where(s.isin(["not_assessed", "no_model"]), "no coverage", s)
    r["reading"] = [READING[(a, b)] for a, b in zip(simp(r.pharm_state), simp(r.tox_state))]
    # direction: only Tox21 distinguishes mode; ChEMBL does not carry agonist/antagonist
    r["direction"] = np.where(r.tox_positive_modes == "", "undetermined", r.tox_positive_modes)
    cols = (["name", "receptor", "reading", "direction"]
            + [c for c in r.columns if c.startswith("pharm_")] + [c for c in r.columns if c.startswith("tox_")])
    return r[cols]


# Voting ceiling per receptor and axis: how many units could vote at most.
def coverage(e):
    """Deployable units per receptor and axis: the vote cap, the same for every molecule."""
    # the unit is the model when by_endpoint ran with unit="model"
    unit_key = "endpoint" if "endpoint" in e.columns else "endpoint_base"
    d = e[e.deployable].drop_duplicates(["receptor", "axis", unit_key])
    t = d.pivot_table(index="receptor", columns="axis", values=unit_key, aggfunc="count", fill_value=0)
    for c in ("pharmacological", "screening"):
        if c not in t:
            t[c] = 0
    t["pharmacological_measures"] = (d[d.axis == "pharmacological"].groupby("receptor").measure
                                   .apply(lambda s: ", ".join(sorted(set(s)))))
    t["screening_modes"] = (d[d.axis == "screening"].groupby("receptor").simple_mode
                          .apply(lambda s: ", ".join(sorted(set(s)))))
    return t.fillna("").reset_index()


if __name__ == "__main__":
    # usage: consensus_v2.py <detail.csv> [output_prefix] [model|endpoint]
    det = sys.argv[1]
    pref = sys.argv[2] if len(sys.argv) > 2 else "consensus_v2"
    unit = sys.argv[3] if len(sys.argv) > 3 else "model"
    d = load(det)
    e = by_endpoint(d, unit=unit)
    r = summarize(e)
    c = coverage(e)
    e.to_csv(f"{pref}_endpoints.csv", index=False)
    r.to_csv(f"{pref}_summary.csv", index=False)
    c.to_csv(f"{pref}_coverage.csv", index=False)
    print(f"voting unit: {unit}")
    print(f"{len(d)} predictions -> {len(e)} units -> {len(r)} molecule x receptor pairs")
    print(r.reading.value_counts().to_string())
