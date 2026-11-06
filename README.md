# EDCQSAR Endocrine Disruption Chemical Quantitative Structure-Activity Relationship

Multilayer consensus predictor that estimates the potential of molecules to interact with nuclear receptors involved in endocrine disruption. The project covers the curation of public data (ChEMBL, BindingDB, Tox21), the training and validation of the models, and the use of the final models to predict new molecules.

## Targets covered

| Receptor | Code |
|---|---|
| Androgen receptor | AR |
| Estrogen receptor alpha / beta | ESR1 / ESR2 |
| Thyroid hormone receptor alpha / beta | THRA / THRB (Tox21: TR, no alpha/beta distinction) |
| PPAR gamma | PPARG |
| Pregnane X receptor | PXR (NR1I2) |

The insulin receptor (IR) is not yet covered.

The **PXR** was added after the others, via `pipeline_pxr.py` (ChEMBL) and `pipeline_pxr_tox21.py` (Tox21).
It is not a hormonal receptor like the others: it is the xenosensor that regulates CYP3A4, and it is included because it mediates part of the disruptors' effect via enzymatic induction and shares ligands with ER, AR, and PPARγ. Three particularities of it, all measured:

- **The assay-type suffix does not indicate what the assay measures.** ChEMBL's `assay_type` field poorly
  separates binding from functional across all targets, but in PXR this is extreme: the largest functional
  group (CYP3A4 reporter in HepG2/DPX-2) is labeled `A` (ADMET), and `B` mixes BODIPY FL vindoline
  displacement, which is true binding, with luciferase reporter. Discarding `A`, as the funnel of the
  other layers does, would throw away 36% of the target's records; merging it with `B` would violate the
  rule of never mixing incompatible assays. The solution adopted was the third: **each `assay_type` is its
  own endpoint**, which is the same rule already applied to `B` and `F`, only extended. No PXR endpoint can
  be called "binding" or "functional" based on the suffix.
- **There is only agonist mode in screening.** The Tox21 program has no PXR antagonism counter-assay, so
  PXR's `direction` never comes out as antagonist. The cytotoxicity mask is **not** applied, because the
  project only applies it to antagonist modes; the `cytotoxic` and `cytotox_confounds` columns are
  calculated and recorded anyway, and 774 of the 2,759 actives (28%) coincide with loss of viability.
- **BindingDB does not add a PXR layer.** There are 2,036 records, but 1,530 are reimports from ChEMBL
  and, after the funnel, 215 measurements remain, with only **31 molecules** absent from the ChEMBL
  layer. The literature fraction, the only one the project's measurement showed to improve class balance,
  is 14 records. Below the funnel's minimum of 50, so no dataset was generated.
## Architecture: three independent layers

The layers are never merged, because they measure different things and the assays are not comparable to each other.

| Layer | Folder | Models | Contents |
|---|---|---|---|
| ChEMBL IC50 | `models_final/` | 14 (+4 archived) | pIC50 regression and classification, by target and assay type (binding `B` or functional `F`) |
| ChEMBL Ki, EC50, and Kd | `models_layers/` | 26 (+14 archived) | same targets, each measure in its own model. PXR has 4 models here (`EC50_A`, `EC50_B`) and 4 archived (`IC50_B`, `EC50_F`); it has no usable Ki or Kd measure |
| Tox21 qHTS | `models_tox21/` | 12 | classification by assay (environmental chemistry). There are 13 curated assays and 12 with a model: the TR agonist one has only 27 actives and is left out. The viability counter-assay only applies to antagonist modes, and the autofluorescence one only to BLA assays |

There are **70 trained models**. The **52 deployable** ones are in the three folders above and are the only ones that the
program loads. The **18 failed** ones are in `models_failed/`, **only as documentation**: the
predictor does not read them, so they do not cost execution time. Each one is described in
`models_failed/README_failed_models.csv`.

### Alert rule

For each molecule-receptor pair, only the layers in which the molecule is **inside the applicability
domain** (maximum Tanimoto similarity to the training set greater than or equal to the model threshold)
and whose model is deployable count. Regarding these layers:

| Alert | Criterion |
|---|---|
| high (pharmacological + environmental) | 2 or more active layers, with at least one ChEMBL and one Tox21 |
| moderate (majority of layers) | fraction of active layers of at least 0.5, with 3 or more layers assessed |
| low (isolated evidence) | at least 1 active layer |
| no alert | no active layer |

A pair in which no layer is within the domain appears as `-` in the table and as `n.d.` in the plots.

## Installation

Requires Python 3.11 or 3.12.

```bash
cd ~/Downloads/pipeline_edc_qsar
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

If `python3 -m venv` fails with `No module named 'ensurepip'` (this happens on Ubuntu, where the module
comes in the `python3-venv` package), use `virtualenv`, which brings `pip` built-in and does not need
`sudo`:

```bash
python3 -m pip install --target .tools virtualenv
PYTHONPATH=.tools python3 -m virtualenv -p python3 .venv
.venv/bin/pip install -r requirements.txt
```

After that the `.tools/` folder can be deleted: it only serves to create the environment.

Installation check:

```bash
# should print 2180
.venv/bin/python -c "import features as F; print(2048 + len(F.DESC_NAMES))"

# should print 'loaded 52' (the 18 failed models stay in models_failed/ and do not enter the count)
.venv/bin/python -c "
import glob, joblib
print('loaded', sum(bool(joblib.load(f)) for f in glob.glob('models_*/*.joblib')))"

# reference prediction: the ERalpha alert should come out 'high'
.venv/bin/python edc.py predict --smiles 'CC(C)(c1ccc(O)cc1)c1ccc(O)cc1'
```

The versions of `rdkit` and `scikit-learn` must be those in `requirements.txt`. With a different version of RDKit the number of descriptors changes and the predictor fails with `Feature shape mismatch`. The `requirements.txt` explains the situation with `xgboost`, whose training version (3.4.2) only exists on conda-forge.

### Shortcut

The `edc` executable in the root calls `edc.py` with the folder's environment, from any directory:

```bash
./edc help
~/Downloads/pipeline_edc_qsar/edc predict --file molecules.csv --name-column name
```
## Usage

The entry point is `edc.py`, with three commands:

```bash
python edc.py help                                               # lists the commands
python edc.py explain                                            # how to use and how to interpret
python edc.py predict --smiles "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1"    # one molecule
python edc.py predict --smiles "SMILES_1" "SMILES_2"              # several molecules
python edc.py predict --file molecules.csv --name-column name     # one list
python edc.py predict --file molecules.csv --axes                 # with the reading in two axes
```

Without arguments, `python edc.py` shows the list of commands. With `--output PREFIX` you choose where
to write; without it, the output goes to `results/YYYYMMDD_HHMMSS_predictions/`. The `explain` command also
counts how many models are installed in each layer.

For the archived run, which in addition to the prediction also generates the plots:

```bash
python run_edc.py --input demo_edcs.csv --name-col name --out demo
```

The input CSV needs a `smiles` column (or specify `--smiles-col`). Each run writes to
`results/YYYYMMDD_HHMMSS_<name>/`:

- `<name>_summary.csv`: alert and layer count per molecule and receptor
- `<name>_detail.csv`: result for each model, with domain, probability or predicted value, and conformal interval
- `<name>_rejected.csv`: molecules rejected during standardization, if any
- `<name>_dashboard.png` and `<name>_active_fraction.png`: alert plots
- `input.csv` and `README.txt`: copy of the input, command, date, and library versions

It is also possible to use the scripts separately:

```bash
python predict_edc2.py --input demo_edcs.csv --name-col name --out demo   # prediction only
python plot_edc_results.py --out demo                                     # plots only, from demo_*.csv
python predict_edc2.py --smiles "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1"           # one molecule
```

`panel_ref.csv` is a reference panel with known classes, and `demo_edcs.csv` is a short example.
## Contents of the repository

| Path | Function |
|---|---|
| `run_edc.py` | runs the prediction, generates the plots, and organizes the output in `results/` |
| `predict_edc2.py` | multilayer consensus predictor |
| `plot_edc_results.py` | alert and active-fraction plots |
| `features.py`, `curation.py` | featurization (ECFP4 by count and RDKit descriptors) and structure standardization |
| `extract_chembl.py`, `extract_bindingdb.py` | extraction of activities from the ChEMBL 37 and BindingDB dumps |
| `pipeline_curate.py`, `pipeline_layers.py`, `pipeline_tox21.py` | curation of the IC50, Ki/Kd/EC50, and Tox21 datasets |
| `modeling.py`, `modeling_layers.py`, `modeling_tox21.py` | training and validation by layer |
| `finalize_models.py`, `conformal.py` | final models and conformal prediction with optimized decision threshold |
| `datasets_layers_by_endpoint/`, `datasets_tox21_by_assay/` | curated Ki/EC50/Kd and Tox21 datasets used in training |
| `models_final/`, `models_layers/`, `models_tox21/` | the 52 deployable models, the only ones the program loads (about 110 MB) |
| `models_failed/` | the 18 failed models, for documentation purposes only, with `README_failed_models.csv` (about 11 MB) |
| `archive_failed_models.py` | moves the failed models to `models_failed/` and returns to execution any that become approved |

Not included: the raw ChEMBL and BindingDB dumps, needed to redo the extraction, and the
curated datasets of the IC50 layer.
## Validation protocol

- External split by Bemis-Murcko scaffold (80/20) and model selection by 5-fold cross-validation on the training set. **Note:** the code fills the *test* set with the **largest** scaffold groups (the docstring of `features.scaffold_split` states the opposite). The external test is, therefore, extrapolation to the dominant chemical families, concentrated in few scaffolds (from 2 to 15 in the Ki/Kd/EC50 endpoints), and the test metrics are a high-variance lower bound. The training cross-validation, in contrast, is random and non-nested, so `Q2_cv` is optimistic.
- Final model refit on 100% of the curated data, with metrics coming from the external scaffold-based test. The metrics recorded inside the production artifact (`Q2_cv`, `R2_test_scaffold`) are, therefore, from **another** fit: the 80% one. The domain threshold recorded in the artifact is also recalculated over 100% of the data, so it is wider than the validated threshold.
- Y-scrambling and applicability domain by Tanimoto similarity.
- Inductive conformal prediction with nominal coverage of 80%, and decision threshold optimized by MCC instead of 0.5. Two measured caveats: the **empirical** coverage falls below 0.80 in 8 of the 9 IC50 endpoints; and the conformal quantile is calibrated on a model fitted on 75% of the training set, but applied to the production model fitted on 100%, such that the coverage guarantee does not hold for the deployed model.
- Curation: valid structures, standardized SMILES, removal of salts and solvents, and deduplication by InChIKey, with the median of replicates and discarding of the compound when they vary by more than 1 log unit. Measure types (IC50, Ki, Kd, EC50) and assay types are never mixed. Stereochemistry is only **recorded** (column `undefined_stereo`, which counts unassigned centers) and never filters. Deduplication uses the full InChIKey, which includes the stereochemical layer, so defined stereoisomers remain as distinct compounds, and a compound with undefined stereochemistry coexists as a separate row from its resolved enantiomers.

## Known issues

Verified in the code and in the data. The breakdown, with the numbers and the diagnostic files, is in
`documents/code_audit_verification.md`.
### Fixed

No fix alters any model's prediction: what changes is which models vote, how
the votes are aggregated, and what the output reports.

**Deployment criterion, now auditable** — `deployment_criterion.py`

The `not_deployable` flag recorded in the `.joblib` files of `models_final/` came from a script missing from
the repository, and was not auditable: `ESR2_B` passed with R² of 0.098 and `AR_F` was failed with 0.076,
with no rule written anywhere. The criterion was reconstructed from the metrics that the artifacts
themselves store:

| Condition | Verdict |
|---|---|
| `n_train < 200` | failed, in any layer |
| regression with `Q2_cv < 0.30` | failed |
| regression with `R2_test_scaffold < 0` | failed (worse than predicting the mean) |
| classification with `ROC_AUC_test < 0.60` | failed |

It reproduces all the recorded flags except two, diverging only on `AR_F` (n = 618, Q² = 0.450,
R² = 0.076, ROC-AUC = 0.772), which becomes deployed. Result: **52 deployed out of 70 trained**,
already including the 9 PXR models (5 deployable: `PXR_EC50_A` and `PXR_EC50_B` in both tasks, plus the
Tox21 assay; `PXR_IC50_B` with n = 176 and `PXR_EC50_F` with n = 62 are failed due to size).
`deployment_criterion.csv` gives the reason for each failure and, in the `local` column, where the model is.

**Failed models are not loaded.** A failed model never votes, so loading it and applying it to every
molecule just to discard the result was cost without use (324 of the panel's 1,260 predictions). They stay in
`models_failed/`, and `predict_edc2.py` ignores with a warning any failed model that reappears in the
execution folders. "Failed" almost always means **little data**, not a bad model: 15 of 18 fail only on
`n < 200`, a project cutoff that was not validated. After retraining, the flow is
`python deployment_criterion.py` and `python archive_failed_models.py --apply`.

Important caveat: **the 12 Tox21 models count as deployed due to absence of a criterion**, not
because they passed one — neither the `.joblib` files nor `README_models_tox21.csv` store a validation
metric. This needs to be stated in any report that uses the screening layer.

**New tier rule, based on units of evidence**

A model is not a piece of evidence. Two redundancies were measured: regression and classification of the same
endpoint use the same compounds and the same features, differing only by the binarization at pAct ≥ 6; and
the Tox21 assays of the same receptor and mode are the same library (up to 99.8% overlap) in different
cell lines. Votes are now collapsed before any counting:

- **pharmacological axis:** 1 unit per ChEMBL endpoint (receptor × measure × assay type). The
  classification decides positivity; the regression enters as potency, not as a separate vote.
- **screening axis:** 1 unit per receptor and mode, collapsing cell lines.
The ladder is now:

The ladder became:

| Level | Condition |
|---|---|
| `not assessed` | no unit assessed, or only one and negative |
| `no alert` | all units assessed negative (2 or more) |
| `high` | the **both** axes positive (convergence) |
| `moderate` | ≥ 2 independent units positive on a single axis |
| `low` | exactly 1 positive unit |

The `not assessed` state is the main improvement: the previous rule called both the molecule
tested and negative and the one that no model managed to assess "no alert". In `panel_ref.csv`, **14 of the 40
pairs that were "no alert" were, in fact, not assessed**. The 10 "high" remained "high", all
by convergence.

**Output**

- The `_detail.csv` brings **only the models that voted**. `--detail deployed` adds those that
  did not vote because they were outside the domain. `--detail all` is accepted and gives the same result, because
  failed models are no longer loaded. The complete table continues to be used internally by
  `--axes`, whose coverage depends on comparing the deployable units with the assessed ones.
- New column `reason_not_voted`: `model failed (<reason>)` or `molecule outside the applicability domain`.
- Every molecule × receptor combination gets a row in `_summary.csv`, including those without any
  coverage, which come out as `not assessed`. Previously, a receptor whose units were all outside the domain
  disappeared from the output.
- The `evidence` column speaks only about what was assessed (`positive 1/3`). The coverage remains in the CSV, in
  the `units_deployable` and `coverage` columns, but is left out of the printed table — the decision was to keep the
  screen showing what was actually measured. It is worth consulting them before comparing two molecules: in the panel, 33
  of the 83 assessed pairs belong to groups of the same level and same positive fraction whose coverage varies
  by more than 0.40 (`Ibuprofen · ERβ` has coverage 0.25 and `17β-Estradiol · TR` has 1.00, and both come out
  "positives 1/1"). The line printed above the table gives the aggregate coverage of the run.
- The `_summary.csv` came to have `evidence`, `units_positive`, `units_assessed`,
  `units_deployable`, `positive_fraction`, `coverage`, `agreement`, `direction` and `alert`.
- The screen reports how many units of evidence were able to vote, out of the total available.

**Produced wrong results without warning**

- **Layer identified by folder name.** The comparison now uses the *base* name
  (`os.path.basename`), so `--model-dirs` with an absolute path works, and a folder with an
  unknown name raises an error. Previously, `--model-dirs /abs/.../models_tox21` labeled the 11 Tox21 models
  as `ChEMBL AR/ESR1/PPARG/TR`, with the assay mode becoming the receptor name.
- **Repeated names in the input.** `predict_edc2.py` now raises an explicit error, instead of merging
  different molecules into one row of the summary.
- **Molecules without a Bemis-Murcko scaffold.** Each one now receives its own group. The external
  test of `TOX21_AR_antagonist_HEK293T` dropped from 95.2% to 0.0% acyclic and that of
  `TOX21_TR_antagonist_GH3` from 91.9% to 0.1%; `PPARG_Kd_B`, whose test had a **single class**,
  came to have two. Affects only retraining.
- **Silent fallback to random split.** `features.scaffold_split` raises an error, and only accepts
  random split with `allow_random=True`.
**They prevented running end to end**

- **BindingDB.** `pipeline_curate.py` reads the `.parquet` that `extract_bindingdb.py` writes and
  **builds** `cache/bindingdb_std_lookup.pkl`, instead of reading a file that no script created.
- **Data folder names.** `modeling_layers.py` and `modeling_tox21.py` accept both `data_layers/`
  and `data_tox21/` as well as `datasets_layers_by_endpoint/` and `datasets_tox21_by_assay/`.
### Open Issues

**Reproducibility of what is deployed** — not fixable by a code change:

- The 18 models in `models_final/` (14 in production and 4 archived in `models_failed/`) **cannot be regenerated by any script in the folder**. The
  `README_calibration.csv` has columns (`production_threshold`, `threshold_gain`, `conformal_mode`,
  `coverage_scaffold`) that `conformal.py` does not produce; it writes another file,
  `results/conformal_and_thresholds.csv`. The script that closes this step is not in the repository, and rewriting it
  would produce equivalent models, not the same ones. (The fifth column, `not_deployable`, no longer depends
  on this script: `deployment_criterion.py` reconstructs it with a written rule.)
- `finalize_models.py` and `conformal.py` depend on `data/dataset_*.csv` (curated IC50 datasets) and on
  `results/model_performance.csv`, both absent. The IC50 datasets were deliberately left out of the
  repository.

**Structural limits of the consensus** — depend on retraining or new data, not on code:

- Votes remain unweighted: two positive units count the same, whether they come from 2,700 compounds
  (`ESR1_B`) or from 232 (`AR_EC50_B`). Weighting would require a confidence measure comparable across the
  layers, which the screening layer does not have.
- The alert cap depends on the receptor: "high" is only attainable for AR, ERα, PPARγ, and PXR, because it requires a
  vote from each layer. ERβ has no Tox21 assay and reaches at most "moderate"; THRβ reaches "low";
  THRα has no deployable ChEMBL model and does not appear in `_summary.csv` (the two-axis reading, with
  `--axes`, shows THRα as `no coverage`, because it knows from the archived models' documentation that it exists). PXR is the thinnest of those that
  reach "high": 2 pharmacological units and 1 screening unit, versus 4 and 4 for ERα. A "no alert" does not mean the same thing across
  different receptors. The `units_deployable` column makes this cap visible.
- **For PXR, the two axes are not on equal footing in the reference panel.** Of the 18 molecules in the
  `panel_ref.csv`, **14 are in the Tox21 PXR training set and none in the ChEMBL PXR training set** (median Tanimoto
  to the nearest training neighbor: 1.00 in screening versus 0.28 in the pharmacological axis). This is why 13 of the 18
  appear as "divergent: screening positive, pharmacological negative" — the screening model recognizes
  compounds it has already seen, while the pharmacological ones are extrapolating. This pattern is a result of the panel's composition, not
  a property of PXR, and the panel is not meant to estimate PXR performance on new molecules.
- "Active" does not mean the same thing in each layer: in ChEMBL it is potency of up to 1 µM with a probability
  threshold of 0.5; in Tox21 it is the depositor's call, with optimized thresholds from 0.11 to 0.40. This is
  why the two axes are never summed into a single fraction.
- The Tox21 layer has no recorded validation metric, so its 12 models have never been submitted to the
  deployment criterion. These are 12 of the 52 units that vote.
- The `panel_ref.csv` **is not a blind evaluation set**: 59.9% of the votes that count come from models
  that already had the molecule in their own training set (84.3% in the screening axis). It serves to check the
  mechanism of the consensus, never to estimate accuracy on new molecules.

The `consensus_v2/` rewrites only the aggregation and makes these points explicit in the report, without retraining
anything. See `documents/consensus_v2_specification.md`.

**Artifact format.** The keys stored inside the `.joblib` model files keep their original names (`limiar`, `modelo`,
`nao_implantavel`, `cobertura_scaffold`, `conformal_modo`, `frac_ativos`, `R2_teste_scaffold`, `ROC_AUC_teste_scaffold`,
`conformal_usa_dificuldade`) and the value `fixo`/`normalizado` of `conformal_modo`, so the 70 existing models stay readable
without re-saving them. Every script that writes or reads these files uses the same names; the columns of the CSV outputs
are in English. Renaming the keys would require re-saving all 70 models (about 127 MB).

## Limitations

- Computational screening. An alert indicates plausibility and does not replace experimental assay or regulatory assessment.
- Molecules outside the applicability domain do not enter the alert rule, even when the model returns a high probability. Compounds from classes poorly represented in training, such as phthalates, may end up without coverage.
- Mechanisms that do not involve the modeled receptors are not detected. For example, effects via aromatase (CYP19A1) or steroidogenesis would require additional endpoints.
- The Tox21 layer is environmental chemistry and uses a different kind of measure (AC50) than the ChEMBL layers.

## Data sources

ChEMBL 37, BindingDB and Tox21. The derived data follow the terms of each source. ChEMBL is distributed under CC BY-SA 3.0. Consult the terms of the other sources before redistributing the datasets.

## License

To be determined.
