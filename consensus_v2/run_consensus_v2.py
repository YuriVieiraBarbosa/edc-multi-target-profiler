#!/usr/bin/env python
"""Shortcut: runs predict_edc2.py and then the two-axis re-aggregation from consensus_v2.py.

Does not modify predict_edc2.py or any model. It just calls both in the correct order, respecting
what predict_edc2.py requires: it identifies the layer by the NAME of the models folder, so it must be
run from the pipeline root and with the default model paths. Passing `--model-dirs` with an
absolute path causes the Tox21 layer to be labeled as ChEMBL and corrupts the votes.

Usage:
    python consensus_v2/run_consensus_v2.py --input molecules.csv --name-col name --out results/run1

Generates:
    <out>_detail.csv, <out>_summary.csv          (from predict_edc2.py, original format)
    <out>_v2_endpoints.csv, <out>_v2_summary.csv,
    <out>_v2_coverage.csv                          (from consensus_v2.py, two axes)
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


# Run predict_edc2.py and then the two-axis re-aggregation, in that order and from the pipeline root.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="Input CSV with the SMILES column")
    ap.add_argument("--smiles-col", default="smiles")
    ap.add_argument("--name-col", default=None)
    ap.add_argument("--out", default="predictions_consensus", help="prefix of the outputs (can be an absolute path)")
    a = ap.parse_args()

    input_data = os.path.abspath(a.input)
    output = os.path.abspath(a.out)
    os.makedirs(os.path.dirname(output) or ".", exist_ok=True)

    cmd = [sys.executable, "predict_edc2.py", "--input", input_data, "--smiles-col", a.smiles_col, "--out", output]
    if a.name_col:
        cmd += ["--name-col", a.name_col]
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    # cwd = pipeline root: the paths models_final/ models_layers/ models_tox21/ are relative and the
    # name of the folder is what defines the layer of each model.
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)

    det = f"{output}_detail.csv"
    if not os.path.exists(det):
        sys.exit(f"predict_edc2.py did not generate {det}")
    subprocess.run([sys.executable, os.path.join(HERE, "consensus_v2.py"), det, f"{output}_v2"], check=True)


if __name__ == "__main__":
    main()
