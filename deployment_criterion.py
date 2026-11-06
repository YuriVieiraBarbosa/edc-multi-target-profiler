#!/usr/bin/env python
"""Explicit deployment criterion for the models.

The `not_deployable` flag recorded inside the .joblib files in models_final/ came from a script that
is not in the repository, so it was not auditable: `ESR2_B` passed with R2 of 0.098 and `AR_F`
was failed with 0.076, with no rule written anywhere. This script reconstructs the
criterion from the metrics that the artifacts themselves store, writes
`deployment_criterion.csv` and becomes the source that predict_edc2.py consults.

Rule, applied per MODEL (not per endpoint):

    n_train < 200                      -> failed, at any layer
    regression:    Q2_cv < 0.30         -> failed
                  R2_test_scaffold<0  -> failed (worse than predicting the mean)
    classification: ROC_AUC_test < 0.60 -> failed

Tox21 is a special case: neither the .joblib files nor the README store a validation metric, so
these 11 models are NOT evaluable by this criterion. They keep voting, but due to
absence of a criterion, not because they passed one - which is recorded in the column
`reason` and must be stated in any report.

Usage: python deployment_criterion.py [output.csv]
"""
import glob
import os
import sys

import joblib
import numpy as np
import pandas as pd

N_MIN = 200
Q2_MIN = 0.30
R2_TEST_MIN = 0.0
ROC_AUC_TEST_MIN = 0.60


FOLDERS = ("models_final", "models_layers", "models_tox21")
ARCHIVE = "models_failed"     # failed models kept only as documentation (see archive_failed_models.py)


# Yield (logical folder, path, where) for every trained .joblib: the run folders and the archive.
# The logical folder stays the ORIGINAL one, because the rule and the filename grammar depend on it.
def _models(root="."):
    for folder in FOLDERS:
        for local, d in (("running", os.path.join(root, folder)),
                         ("archive", os.path.join(root, ARCHIVE, folder))):
            for p in sorted(glob.glob(os.path.join(d, "*.joblib"))):
                yield folder, p, local


# Read every trained model (running or archived) and collect the metrics each layer stores.
def collect(root="."):
    cam = os.path.join(root, "models_layers", "README_models_layers.csv")
    lay = pd.read_csv(cam).set_index("endpoint") if os.path.exists(cam) else pd.DataFrame()
    rows = []
    for folder, p, local in _models(root):
        m = joblib.load(p)
        base = os.path.basename(p)[:-7]
        # the 'task' field is not consistent across layers (there is 'regression' and 'regressao'),
        # so the task comes from the filename suffix, which is uniform across models
        reg = base.endswith("_regression")
        r = dict(folder=folder, model_id=base, local=local,
                 task="regression" if reg else "classification",
                 # the keys of the saved .joblib artifacts keep their original names (limiar, modelo, nao_implantavel, ...) so the existing model files stay readable
                 saved_flag=bool(m.get("nao_implantavel", False)),
                 n=np.nan, cv_metric=np.nan, test_metric=np.nan)
        if folder == "models_final":
            r["n"] = m.get("n_train")
            r["cv_metric"] = m.get("Q2_cv") if reg else m.get("ROC_AUC_cv")
            r["test_metric"] = (m.get("R2_teste_scaffold") if reg
                                  else m.get("ROC_AUC_teste_scaffold"))
        elif folder == "models_layers":
            ep = base.rsplit("_", 1)[0]
            if ep in lay.index:
                r["n"] = lay.loc[ep, "n"]
                if reg:
                    r["cv_metric"] = lay.loc[ep, "Q2_cv"]
        else:
            r["n"] = m.get("n_train")
        rows.append(r)
    return pd.DataFrame(rows)


# Apply the rule above to one model, returning (rejected?, reason).
def judge(r):
    if r.folder == "models_tox21":
        return False, "no validation metric recorded (not assessable)"
    if pd.notna(r.n) and r.n < N_MIN:
        return True, f"n={int(r.n)} < {N_MIN}"
    if r.task == "regression":
        if pd.notna(r.cv_metric) and r.cv_metric < Q2_MIN:
            return True, f"Q2_cv={r.cv_metric:.3f} < {Q2_MIN:.2f}"
        if pd.notna(r.test_metric) and r.test_metric < R2_TEST_MIN:
            return True, f"R2_test_scaffold={r.test_metric:.3f} < {R2_TEST_MIN:.2f}"
    elif pd.notna(r.test_metric) and r.test_metric < ROC_AUC_TEST_MIN:
        return True, f"ROC_AUC_test={r.test_metric:.3f} < {ROC_AUC_TEST_MIN:.2f}"
    return False, "passed"


# Build the full table: one row per model, with the verdict and the reason behind it.
def table(root="."):
    T = collect(root)
    T[["not_deployable", "reason"]] = T.apply(lambda r: pd.Series(judge(r)), axis=1)
    T["diverged_from_saved_flag"] = T.not_deployable != T.saved_flag
    # the order cannot depend on where the file is today, otherwise moving a model
    # rewrites the entire CSV and the git diff stops showing what actually changed
    T = T.sort_values(["folder", "model_id"]).reset_index(drop=True)
    cols = [c for c in T.columns if c != "local"] + ["local"]
    return T[cols]


if __name__ == "__main__":
    output = sys.argv[1] if len(sys.argv) > 1 else "deployment_criterion.csv"
    T = table()
    T.to_csv(output, index=False)
    n_rep = int(T.not_deployable.sum())
    print(f"{len(T)} models | deployed: {len(T) - n_rep} | failed: {n_rep}")
    misplaced = T[(T.local == "running") & T.not_deployable]
    if len(misplaced):
        print(f"  WARNING: {len(misplaced)} failed model(s) still in the running folders; "
              "run archive_failed_models.py --apply")
    print(f"  of deployed, {int((T.folder == 'models_tox21').sum())} are from Tox21 and were not "
          "assessed due to lack of a recorded metric")
    div = T[T.diverged_from_saved_flag]
    if len(div):
        print(f"\ndivergences from the flag stored in the .joblib file ({len(div)}):")
        for _, r in div.iterrows():
            print(f"  {r.model_id:28s} stored={'failed' if r.saved_flag else 'passed':9s}"
                  f" criterion={'failed' if r.not_deployable else 'passed':9s} ({r.reason})")
    print(f"\nwritten to {output}")
