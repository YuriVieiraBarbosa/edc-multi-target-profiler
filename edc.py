#!/usr/bin/env python
"""Entry point for the EDC-QSAR pipeline.

Three commands:
    predict   predicts a molecule or a list of molecules
    explain   how to use the program and how to read the result
    help      lists the commands above

Examples:
    python edc.py predict --smiles "CC(=O)Oc1ccccc1C(=O)O"
    python edc.py predict --file molecules.csv --name-column name --axes
    python edc.py explain
    python edc.py help
"""
import argparse
import datetime
import glob
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_FOLDERS = ("models_final", "models_layers", "models_tox21")
LAYER_LABEL = {"models_final": "ChEMBL IC50", "models_layers": "ChEMBL Ki/Kd/EC50",
                 "models_tox21": "Tox21 (screening)"}


# Count the .joblib files of each model directory, without loading them.
def _inventory():
    """Per folder: (trained, deployed). Trained is not deployed - the failed models
    per deployment_criterion.csv remain in models_failed/, only as documentation, and
    are never loaded. 'Trained' sums both folders; 'deployed' only the running folder."""
    import pandas as pd
    p = os.path.join(HERE, "deployment_criterion.csv")
    rep = set()
    if os.path.exists(p):
        t = pd.read_csv(p)
        rep = set(t.loc[t.not_deployable, "model_id"])
    inv = {}
    for d in MODEL_FOLDERS:
        ids = [os.path.basename(f)[:-7] for f in glob.glob(os.path.join(HERE, d, "*.joblib"))]
        fname = glob.glob(os.path.join(HERE, "models_failed", d, "*.joblib"))
        inv[d] = (len(ids) + len(fname), sum(1 for i in ids if i not in rep))
    return inv


# Receptors that were trained but have no deployable model, read from the archive's documentation.
def _receptors_without_model(det):
    import pandas as pd
    p = os.path.join(HERE, "models_failed", "README_failed_models.csv")
    if not os.path.exists(p):
        return []
    # the receptor is the first segment of the model id, the same rule as predict_edc2.load
    documented = {m.split("_")[0] for m in pd.read_csv(p).model_id}
    return sorted(documented - set(det.receptor))


# Resolve the output prefix; with no --output, create results/<timestamp>_<name>/.
def _output_prefix(output, default_name="predictions"):
    if output:
        destination = os.path.abspath(output)
        os.makedirs(os.path.dirname(destination) or ".", exist_ok=True)
        return destination
    now = datetime.datetime.now()
    folder = os.path.join(os.getcwd(), "results", f"{now:%Y%m%d_%H%M%S}_{default_name}")
    os.makedirs(folder, exist_ok=True)
    return os.path.join(folder, default_name)


# Read the molecules from --smiles or from a CSV file, returning SMILES and names.
def _read_input(a, ap):
    if a.smiles:
        # without --file there is no name column; labeling mol_1..mol_n avoids a table
        # with the full SMILES on each row, and the legend below preserves the link
        names = [f"mol_{i + 1}" for i in range(len(a.smiles))]
        print("Molecules received:")
        for n, s in zip(names, a.smiles):
            print(f"  {n} = {s}")
        return list(a.smiles), names
    if a.file:
        if not os.path.exists(a.file):
            ap.error(f"file not found: {a.file}")
        df = pd.read_csv(a.file)
        if a.smiles_column not in df.columns:
            ap.error(f"column '{a.smiles_column}' does not exist in {a.file}. "
                     f"Available columns: {', '.join(df.columns)}")
        names = None
        if a.name_column:
            if a.name_column not in df.columns:
                ap.error(f"column '{a.name_column}' does not exist in {a.file}")
            names = df[a.name_column].astype(str).tolist()
        return df[a.smiles_column].astype(str).tolist(), names
    ap.error("provide --smiles or --file")


# ANSI colours, by alert level and by two-axis reading.
RESET = "\033[0m"
ALERT_COLOR = {"high": "\033[1;97;41m",      # white on red
              "moderate": "\033[1;33m",      # yellow
              "low": "\033[36m",           # cyan
              "no alert": "\033[90m",      # gray
              "not assessed": "\033[35m"}    # magenta: absence of evidence, not evidence of absence
READING_COLOR = (("convergent positive", "\033[1;97;41m"),
               ("divergent", "\033[1;35m"),
               ("positive", "\033[1;33m"),
               ("no coverage", "\033[90m"),
               ("negative", "\033[90m"))


# Decide whether to emit ANSI colour: never when piped to a file or when the user opted out.
def _active_color(no_color):
    if no_color or os.environ.get("NO_COLOR"):
        return False
    return bool(getattr(sys.stdout, "isatty", lambda: False)())


# Pick the colour of one cell: exact alert level, or first matching substring of a reading.
def _color_for(value, column):
    text = str(value)
    if column == "alert":
        return ALERT_COLOR.get(text.split(" (")[0], "")
    for key, fill_color in READING_COLOR:
        if key in text:
            return fill_color
    return ""


# Render a table, colouring one column. Widths come from the PLAIN text, and the escape
# codes are added after padding, so the alignment does not count invisible characters.
def _show(df, color_column=None, use_color=False, limit=25, file=""):
    vis = df.head(limit)
    cols = list(vis.columns)
    # float with 2 decimal places; the remainder as text. str(1.0) would output "1.0" and would break the column.
    def _fmt(v):
        if isinstance(v, float):
            # v != v is NaN: appears when nothing was assessed (0 of 0), and "nan" on the screen
            # looks like a defect instead of absence
            return f"{v:.2f}" if v == v else "-"
        return str(v)
    txt = {c: [_fmt(v) for v in vis[c]] for c in cols}
    width = {c: max([len(c)] + [len(v) for v in txt[c]]) for c in cols}
    print("  ".join(c.rjust(width[c]) for c in cols))
    for i in range(len(vis)):
        cells = []
        for c in cols:
            raw_records = txt[c][i]
            fill_color = _color_for(raw_records, c) if (use_color and c == color_column) else ""
            # the fill stays OUTSIDE the color code, otherwise the colored background of the
            # level "high" extends through the spaces and becomes a bar.
            cells.append(" " * (width[c] - len(raw_records)) + f"{fill_color}{raw_records}{RESET}" if fill_color
                           else raw_records.rjust(width[c]))
        print("  ".join(cells))
    if len(df) > limit:
        print(f"... and {len(df) - limit} more row(s). The complete table is at {file}")


# One-line colour legend, printed only when colour is on.
def _legend_colors(column):
    if column == "alert":
        items = [(n, ALERT_COLOR[n]) for n in ("high", "moderate", "low", "no alert", "not assessed")]
    else:
        items = [("convergent positive", READING_COLOR[0][1]), ("divergent", READING_COLOR[1][1]),
                 ("only one positive axis", READING_COLOR[2][1]), ("negative / no coverage", READING_COLOR[3][1])]
    print("legend: " + "   ".join(f"{fill_color}{name}{RESET}" for name, fill_color in items))


# COMMAND 1: predict one molecule or a list of molecules and print the consensus table.
def predict_cmd(a, ap):
    fill_color = _active_color(a.no_color)
    smi, names = _read_input(a, ap)
    # predict_edc2.py identifies the layer by the NAME of the models folder, with a relative path;
    # that's why the prediction runs from the root of the pipeline.
    previous = os.getcwd()
    sys.path.insert(0, HERE)
    os.chdir(HERE)
    try:
        import predict_edc2
        det, res = predict_edc2.predict(smi, names)
    finally:
        os.chdir(previous)

    prefix = _output_prefix(a.output)
    det_complete = det
    det = predict_edc2.filter_detail(det, a.detail)
    det.to_csv(f"{prefix}_detail.csv", index=False)
    res.to_csv(f"{prefix}_summary.csv", index=False)

    # Molecules that were not scored: the structural standardization rejected them
    # (invalid SMILES, inorganic element, salt without an organic fragment, mass outside
    # 100-1000 Da). predict_edc2.predict only returns the reason when ALL fail, so here
    # the reason is recalculated only for the ones missing -- cheap, and makes the message useful.
    input_names = list(names) if names is not None else list(smi)
    scored = set(det.name.unique()) if "endpoint" in det.columns else set()
    missing = [(n, s) for n, s in zip(input_names, smi) if n not in scored]
    rej = pd.DataFrame(columns=["name", "input_smiles", "rejection_reason"])
    if missing:
        sys.path.insert(0, HERE)
        import curation
        rej = pd.DataFrame([{"name": n, "input_smiles": s,
                             "rejection_reason": curation.standardize(s)["reject_reason"]}
                            for n, s in missing])
        rej.to_csv(f"{prefix}_rejected.csv", index=False)

    print(f"\n{len(smi)} molecule(s) submitted; {len(scored)} scored by the models.")
    if len(rej):
        print(f"{len(rej)} rejected during structural standardization:")
        for r in rej.itertuples():
            print(f"   {r.name}: {r.rejection_reason}")
        print(f"   (detail in {prefix}_rejected.csv)")
    if len(scored) == 0:
        print("No molecule reached the models. Fix the input and run again.")
    elif len(res) == 0:
        print("No model can vote: all the molecules fell outside the domain "
              "of applicability of all deployable models.")
    else:
        columns = ["name", "receptor", "evidence", "positive_fraction",
                   "agreement", "direction", "alert"]
        print()
        n_voted_units = int(res.units_assessed.sum())
        cap = int(res.units_deployable.sum())
        is_deployable = int((~det_complete["not_deployable"]).sum())
        n_voting_models = int(det_complete["counted"].sum())
        print(f"{n_voting_models} of {is_deployable} deployed models voted, forming {n_voted_units} of {cap} "
              f"evidence units ({n_voted_units / cap:.0%} of coverage).")
        if a.detail == "voters":
            print(f"The _detail.csv brings only these {n_voting_models}. Use --detail deployed "
                  "to see those that were outside the domain and why.")
        _show(res[columns], color_column="alert", use_color=fill_color,
                 file=f"{prefix}_summary.csv")
        if fill_color:
            _legend_colors("alert")
        without_coverage = res[res.units_assessed == 0]
        if len(without_coverage):
            print(f"\nWarning: {len(without_coverage)} molecule-receptor pair(s) with no vote "
                  "(molecule outside the domain of all models for that receptor).")

    if a.axes:
        path_v2 = os.path.join(HERE, "consensus_v2", "consensus_v2.py")
        if not os.path.exists(path_v2):
            print("\n[--axes ignored] consensus_v2/consensus_v2.py not found.")
        else:
            sys.path.insert(0, os.path.dirname(path_v2))
            import consensus_v2
            # the FULL table, because v2's coverage compares deployable against assessed
            e = consensus_v2.by_endpoint(consensus_v2.load(det_complete), unit=a.unit)
            r = consensus_v2.summarize(e, receptors_without_model=_receptors_without_model(det_complete))
            r.to_csv(f"{prefix}_axes_summary.csv", index=False)
            e.to_csv(f"{prefix}_axes_units.csv", index=False)
            print(f"\nReading on two axes (vote unit: {a.unit}):")
            _show(r[["name", "receptor", "reading", "direction"]], color_column="reading",
                     use_color=fill_color, file=f"{prefix}_axes_summary.csv")
            if fill_color:
                _legend_colors("reading")

    print(f"\nFiles saved with prefix: {prefix}")
    return 0


EXPLANATION_TEXT = """
==============================================================================
 EDC-QSAR - how to use it and how to read the result
==============================================================================

WHAT THE PROGRAM DOES
  It takes a molecule (SMILES) and returns, for each receptor, an indication of
  possible endocrine disruption. It does NOT measure anything: it is a
  computational screen for prioritising compounds, never proof of activity.

  Receptors covered: AR, ERalpha (ESR1), ERbeta (ESR2), PPARgamma, THRalpha,
  THRbeta and PXR (pregnane X receptor, NR1I2).

  PXR is not a hormone receptor like the others: it is the xenosensor that
  regulates CYP3A4. It is included because it mediates part of the effect of
  endocrine disruptors through enzyme induction and shares ligands with ER, AR
  and PPARgamma. For PXR the screening layer only has an agonist mode, so its
  direction never comes out as antagonist, and the assay-type suffix (A, B, F)
  does not tell you whether the assay is a binding or a functional one.

HOW TO USE IT
  1) One molecule:
       python edc.py predict --smiles "CC(=O)Oc1ccccc1C(=O)O"
  2) Several molecules at once:
       python edc.py predict --smiles "SMILES_1" "SMILES_2" "SMILES_3"
  3) A CSV file (recommended for lists):
       python edc.py predict --file molecules.csv --name-column name
     The CSV needs a SMILES column (default: 'smiles'). Use
     --smiles-column if it has another name, and --name-column for the label.
  4) With the two-axis reading, which keeps pharmacological evidence apart from
     environmental screening:
       python edc.py predict --file molecules.csv --axes

HOW THE RESULT IS FORMED
  Every deployable model casts one binary vote (active / inactive):
    - regression      votes active if the predicted value is >= 6.0 (potency <= 1 uM)
    - classification  votes active if the probability passes the model's threshold
  Only a model that is deployable AND whose applicability domain contains the
  molecule votes. A model outside its domain does NOT vote, and it also does NOT
  count as an inactive vote - which is why coverage is shown next to the level.

  By default the _detail.csv contains ONLY the models that voted. To also see the
  ones that fell outside the domain, and why (column reason_not_voted), use
  --detail deployed. Failed models are not loaded: they stay in models_failed/
  as documentation only, with the reason for each one.

UNITS OF EVIDENCE - why a model is not a piece of evidence
  Counting one vote per model overstates the evidence, because of two
  redundancies measured in this project:
    - regression and classification of the same endpoint use the SAME compounds
      and the same features, differing only in the binarisation at pAct >= 6;
    - the Tox21 assays for one receptor and mode are the SAME compound library
      (up to 99.8% in common) run in different cell lines.
  So the votes are first collapsed into UNITS:
    - pharmacological axis: 1 unit per ChEMBL endpoint (receptor x measure x
      assay type). The classification decides whether the unit is positive; the
      regression enters as the predicted potency, not as a separate vote.
    - screening axis: 1 unit per receptor and mode (agonist / antagonist),
      collapsing the cell lines. The unit is positive if any cell line flags it.

THE TWO AXES
  The axis says WHERE the evidence comes from: ChEMBL (pharmacological: affinity
  and potency in a target-directed assay) or Tox21 (screening: cell-based qHTS).
  The two do not substitute for each other, which is why the highest level
  requires both.

COLUMNS OF _summary.csv
  evidence              e.g. "positive 1/3": 1 positive unit out of the 3
                        that could be assessed
  units_positive        positive units of evidence
  units_assessed        units that could be assessed
  units_deployable      the receptor's ceiling: units that exist
  positive_fraction     positive / assessed
  coverage              assessed / deployable  (0 to 1)
  pharm_units_*         the same, on the pharmacological axis only
  screen_units_*        the same, on the screening axis only
  agreement             convergent, divergent, single axis or negative
  direction             agonist / antagonist, from the ACTIVE Tox21 votes
  alert                 high, moderate, low, no alert or not assessed

  The printed table shows only what was assessed. The columns 'coverage' and
  'units_deployable' are not shown on screen but are in the _summary.csv, and
  it is worth checking them before comparing two molecules: a positive
  fraction of 1.00 with coverage 0.17 means that a single unit out of seven
  gave an opinion, not that the evidence is unanimous. The line printed just
  above the table gives the aggregate coverage of the run.

ALERT LEVELS (the first rule that applies wins)
  not assessed  no unit assessed, or a single unit that is negative -
                insufficient coverage to claim anything
  no alert      all assessed units are negative (2 or more)
  high          BOTH axes positive (convergence)
  moderate      2 or more independent positive units on a single axis
  low           exactly 1 positive unit

  'not assessed' is a state that the previous rule did not have: it used the
  label "no alert" both for a molecule that was tested and came out negative
  and for one that no model could assess. Those are opposite things.

HOW TO INTERPRET IT - what the alert does NOT say
  1. The direction comes ONLY from Tox21, which is the layer that records the
     assay mode. On the pharmacological axis it stays "undetermined", and the
     level itself does not separate agonism from antagonism.
  2. The ceiling depends on the receptor. "High" requires convergence between
     the two axes, so it is only reachable where both exist: AR, ERalpha,
     PPARgamma and PXR - the last one being the thinnest of the four, with 2
     pharmacological units and 1 screening unit, against 4 and 4 for ERalpha.
     ERbeta has no Tox21 assay and stops at moderate; THRalpha has no deployable
     model and does not appear. The units_deployable column shows this ceiling,
     which is why it should be read together with the level.
  3. "Active" does not mean the same thing in each layer. In ChEMBL it is a
     potency of up to 1 uM; in Tox21 it is the assay's active/inactive call,
     with a more permissive probability threshold.
  4. The votes are not independent. Ki and EC50 for the same receptor share
     compounds, and the Tox21 assays for a receptor are the same modes in
     different cell lines.
  5. Close to the 6.0 cutoff the regression vote is noise. Prefer the
     classification probability as the indicator of activity, and read the
     regression as an order of magnitude of potency, together with the
     conformal interval (columns interval_lower and interval_upper of the
     _detail.csv).
  6. Low coverage is common for drug-like molecules: the ChEMBL models are
     narrow medicinal-chemistry series, so many drugs fall outside their
     domain and the verdict ends up coming from Tox21 alone.

FILES WRITTEN
  <prefix>_summary.csv       one row per molecule x receptor
  <prefix>_detail.csv        one row per model, with probability,
                             predicted value, conformal interval, similarity
                             to the training set and whether the vote counted
  with --axes:
  <prefix>_axes_summary.csv  two-axis reading, without summing the layers
  <prefix>_axes_units.csv    the evidence per voting unit

OTHER PROGRAMS IN THE FOLDER
  python run_edc.py --input list.csv --name-col name --out run
      same prediction, plus the charts, archived in results/
  python consensus_v2/consensus_v2.py <detail.csv> <prefix> [model|endpoint]
      only the two-axis re-aggregation, from an existing _detail.csv
==============================================================================
"""


# COMMAND 2: print the usage and interpretation guide, plus a live model inventory.
def explain_cmd(a, ap):
    print(EXPLANATION_TEXT.strip())
    inv = _inventory()
    total = sum(t for t, _ in inv.values())
    is_deployable = sum(i for _, i in inv.values())
    print("\nMODELS INSTALLED IN THIS FOLDER  (deployed / trained)")
    for d, (t, i) in inv.items():
        print(f"  {LAYER_LABEL[d]:22s} {i:3d} / {t:<3d}  ({d}/)")
    print(f"  {'TOTAL':22s} {is_deployable:3d} / {total:<3d}")
    print("  Trained but not deployed: failed ones remain in models_failed/ only as")
    print("  documentation and are never loaded.")
    print("  Tox21 ones count as deployed due to ABSENCE of validation metric")
    print("  recorded, not because they passed a criterion. Reason per model in")
    print("  criterion_deployment.csv.")
    if total == 0:
        print("\n  WARNING: no model found. Run this command from the root of the"
              "\n  pipeline, or check if the models_* folders are in place.")
    return 0


# COMMAND 3: list the available commands and the most common options.
def help_cmd(a, ap):
    print("""
EDC-QSAR - available commands

  1) predict     predicts a molecule or a list of molecules
       python edc.py predict --smiles "C1=CC=CC=C1O"
       python edc.py predict --smiles "SMILES_1" "SMILES_2"
       python edc.py predict --file molecules.csv --name-column name
       python edc.py predict --file molecules.csv --axes --output output/test1

     options:
       --smiles ...           one or more SMILES on the command line
       --file CSV          file with a SMILES column
       --smiles-column NAME   name of the SMILES column   (default: smiles)
       --name-column NAME     column with the molecule label  (optional)
       --output PREFIX        prefix of the output files
                              (default: results/<timestamp>_predictions)
       --axes                also records the reading on two axes
       --unit model|endpoint   voting unit used with --axes (default: model)
       --no-color              turns off the colors (already off when redirecting)
       --detail MODE         rows of _detail.csv: voters (default),
                              deployed or all

  2) explain   how to use the program and how to read the result
       python edc.py explain

  3) help      shows this list
       python edc.py help

Detail of each command:  python edc.py <command> --help
""".strip())
    return 0


# Build the argument parser with one subcommand per function.
def build_parser():
    ap = argparse.ArgumentParser(
        prog="edc.py", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", metavar="command")

    p1 = sub.add_parser("predict", help="predicts a molecule or a list of molecules",
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    p1.add_argument("--smiles", nargs="+", help="one or more SMILES")
    p1.add_argument("--file", help="CSV with a SMILES column")
    p1.add_argument("--smiles-column", default="smiles", dest="smiles_column",
                    help="name of the SMILES column in the CSV (default: smiles)")
    p1.add_argument("--name-column", default=None, dest="name_column",
                    help="column with the molecule label (optional)")
    p1.add_argument("--output", default=None,
                    help="prefix of the output files (default: results/<timestamp>_predictions)")
    p1.add_argument("--axes", action="store_true",
                    help="also records the reading on two axes (pharmacological x screening)")
    p1.add_argument("--unit", choices=("model", "endpoint"), default="model",
                    help="voting unit used with --axes (default: model)")
    p1.add_argument("--detail", choices=("voters", "deployed", "all"),
                    default="voters",
                    help="which models go to _detail.csv (default: voters)")
    p1.add_argument("--no-color", action="store_true", dest="no_color",
                    help="turns off the terminal colors (already off when redirecting the output)")
    p1.set_defaults(func=predict_cmd)

    p2 = sub.add_parser("explain", help="how to use the program and how to read the result")
    p2.set_defaults(func=explain_cmd)

    p3 = sub.add_parser("help", help="lists the available commands")
    p3.set_defaults(func=help_cmd)
    return ap


# Entry point: with no command, show the command list instead of an error.
def main(argv=None):
    ap = build_parser()
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        return help_cmd(None, ap)
    a = ap.parse_args(argv)
    if not getattr(a, "func", None):
        return help_cmd(a, ap)
    return a.func(a, ap)


if __name__ == "__main__":
    sys.exit(main())
