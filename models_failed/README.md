# Failed models (documentation only)

These models were trained and **are not loaded by the predictor**: `predict_edc2.py`, `edc.py`
and `run_edc.py` only look at `models_final/`, `models_layers/` and `models_tox21/`. They are here to
record what was attempted and why it was not included. Each one has the corresponding row in
`README_failed_models.csv` and in `../deployment_criterion.csv`.

"Failed" almost always means **little data**, not a bad model: the criterion requires n >= 200 and
most of the failed models fail only on that. The cutoff of 200 is a project convention, not a
validated value. See the `reason` column.

The indexes `models_final/README_calibration.csv` and `models_layers/README_models_layers.csv`
continue to list all trained endpoints, including these, so that an endpoint listed
there may not have a `.joblib` file in the running folder.

## How a model returns to running

Nothing is moved by hand. Retrain with more data, run `python deployment_criterion.py` and then
`python archive_failed_models.py --apply`: whatever passes the criterion returns to the origin folder and
whatever fails comes here. The `.joblib` file keeps the same name in both places.
