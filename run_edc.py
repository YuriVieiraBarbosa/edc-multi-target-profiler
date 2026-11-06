#!/usr/bin/env python
"""Runs predict_edc2.py (no changes), generates the plots and stores everything in
its own folder inside results/.

Accepts the same arguments as predict_edc2.py:
    python run_edc.py --input demo_edcs.csv --name-col name --out demo

Each running creates:  results/YYYYMMDD_HHMMSS_<out>/
    <out>_summary.csv, <out>_detail.csv, <out>_rejected.csv (if any),
    <out>_dashboard.png, <out>_active_fraction.png,
    input.csv (copy of the input file) and README.txt (command, date, versions).
"""
import argparse
import datetime
import os
import platform
import shutil
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(add_help=False)
ap.add_argument("--out", default="predictions_consensus")
ap.add_argument("--input")
a, remainder = ap.parse_known_args()

if "-h" in remainder or "--help" in remainder:
    sys.exit(subprocess.call([sys.executable, os.path.join(here, "predict_edc2.py")] + remainder))

name = os.path.basename(os.path.normpath(a.out)) or "predictions_consensus"
now = datetime.datetime.now()
result_base = os.path.join(os.getcwd(), "results")
os.makedirs(result_base, exist_ok=True)
folder, k = None, 0
while folder is None:
    cand = os.path.join(result_base, f"{now:%Y%m%d_%H%M%S}_{name}" + (f"_{k}" if k else ""))
    try:
        os.makedirs(cand)
        folder = cand
    except FileExistsError:
        k += 1
prefix = os.path.join(folder, name)

cmd = [sys.executable, os.path.join(here, "predict_edc2.py"), "--out", prefix]
if a.input:
    cmd += ["--input", a.input]
cmd += remainder
rc = subprocess.call(cmd)
if rc != 0:
    if not os.listdir(folder):
        os.rmdir(folder)
    sys.exit(rc)

if a.input and os.path.exists(a.input):
    shutil.copy2(a.input, os.path.join(folder, "input.csv"))

sys.path.insert(0, here)
import plot_edc_results
plot_edc_results.generate_all(prefix)


# Installed version of a module, or '?' when it cannot be imported.
def _version(mod):
    try:
        return __import__(mod).__version__
    except Exception:
        return "?"


rows = [
    f"Running: {now:%Y-%m-%d %H:%M:%S}",
    "Command: python run_edc.py " + " ".join(sys.argv[1:]),
    f"Working folder: {os.getcwd()}",
    f"Python {platform.python_version()} | rdkit {_version('rdkit')} | "
    f"scikit-learn {_version('sklearn')} | pandas {_version('pandas')} | numpy {_version('numpy')}",
    "",
    "Files:",
] + [f"  {f}" for f in sorted(os.listdir(folder))]
with open(os.path.join(folder, "README.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(rows) + "\n")

print(f"\nResults of this running: {os.path.relpath(folder)}")
for f in sorted(os.listdir(folder)):
    print(f"   - {f}")
