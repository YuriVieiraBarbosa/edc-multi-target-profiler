#!/usr/bin/env python
"""Keeps failed models out of running, stored only as documentation.

A model failed by deployment_criterion.csv never votes. Leaving it in models_final/,
models_layers/ or models_tox21/ only costs time: the predictor loads it and applies it to every
molecule only to then discard the result. This script moves the failed models to
models_failed/<origin folder>/, from where the predictor loads nothing, and keeps the
record of each one in models_failed/README_failed_models.csv.

It is symmetric: an archived model that becomes APPROVED (for example after
retraining with more data and running deployment_criterion.py again) goes back to the running
folder. No model is deleted or altered, only moved between folders.

Workflow after training new models:
    python deployment_criterion.py          # rewrites the criterion, looking at both folders
    python archive_failed_models.py           # shows what would move between folders
    python archive_failed_models.py --apply # actually moves them

Usage: python archive_failed_models.py [--apply]
"""
import os
import shutil
import sys

import pandas as pd

ARCHIVE = "models_failed"
CRITERION = "deployment_criterion.csv"

README_TEXT = """# Failed models (documentation only)

These models were trained and **are not loaded by the predictor**: `predict_edc2.py`, `edc.py`
and `run_edc.py` only look at `models_final/`, `models_layers/`, and `models_tox21/`. They are here to
record what was tried and why it did not make it in. Each one has the corresponding row in
`README_failed_models.csv` and in `../deployment_criterion.csv`.

"Failed" almost always means **too little data**, not a bad model: the criterion requires n >= 200 and
most of the failed models fail only on that. The cutoff of 200 is a project convention, not a
validated value. See the `reason` column.

The indexes `models_final/README_calibration.csv` and `models_layers/README_models_layers.csv`
keep listing all the trained endpoints, including these, so an endpoint listed
there may not have a `.joblib` in the running folder.

## How a model goes back to running

There is no manual moving. Retrain with more data, run `python deployment_criterion.py` and then
`python archive_failed_models.py --apply`: whatever passes the criterion goes back to the origin folder and
whatever fails comes here. The `.joblib` file keeps the same name in both places.
"""


# Compare where each model is with where the criterion says it belongs; return the moves.
def plan(root="."):
    t = pd.read_csv(os.path.join(root, CRITERION))
    move = []
    for r in t.itertuples():
        exe = os.path.join(root, r.folder, f"{r.model_id}.joblib")
        fname = os.path.join(root, ARCHIVE, r.folder, f"{r.model_id}.joblib")
        if r.not_deployable and os.path.exists(exe):
            move.append((r.model_id, "archive", exe, fname, r.reason))
        elif not r.not_deployable and os.path.exists(fname):
            move.append((r.model_id, "promote", fname, exe, r.reason))
    return t, move


# Write the documentation of the archived models, from the criterion table.
def write_docs(t, root="."):
    rep = t[t.not_deployable][["model_id", "folder", "task", "n", "cv_metric",
                                "test_metric", "reason"]]
    rep = rep.rename(columns={"folder": "origin_folder"})
    os.makedirs(os.path.join(root, ARCHIVE), exist_ok=True)
    rep.to_csv(os.path.join(root, ARCHIVE, "README_failed_models.csv"), index=False)
    with open(os.path.join(root, ARCHIVE, "README.md"), "w") as fh:
        fh.write(README_TEXT)
    return len(rep)


if __name__ == "__main__":
    apply = "--apply" in sys.argv[1:]
    if not os.path.exists(CRITERION):
        raise SystemExit(f"{CRITERION} not found. Run python deployment_criterion.py first.")
    t, move = plan()
    n_files = sum(1 for m in move if m[1] == "archive")
    n_pro = len(move) - n_files
    print(f"{len(t)} models in the criterion | to archive: {n_files} | to promote: {n_pro}")
    conflict = [m for m in move if os.path.exists(m[3])]
    if conflict:
        for m in conflict:
            print(f"  CONFLICT: {m[0]} exists in both places; nothing was moved.")
        raise SystemExit(1)
    for name, action, origin, destination, reason in move:
        print(f"  {action:8s} {name:28s} ({reason})")
    if not apply:
        print("\nsimulation: nothing was moved. Use --apply to move.")
    else:
        for name, action, origin, destination, reason in move:
            os.makedirs(os.path.dirname(destination), exist_ok=True)
            shutil.move(origin, destination)
        n = write_docs(t)
        print(f"\nmoved {len(move)} | {n} models documented in {ARCHIVE}/README_failed_models.csv")
