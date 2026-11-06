# QSAR Pipeline for Endocrine Disruptors

Multilayer consensus predictor that estimates the potential of molecules to interact with nuclear receptors
involved in endocrine disruption. The project covers the curation of public data
(ChEMBL, BindingDB, Tox21), the training and validation of the models, and the use of the final models to
predict new molecules.

## Targets Covered

| Receptor | Code |
|---|---|
| Androgen receptor | AR |
| Estrogen receptor alpha / beta | ESR1 / ESR2 |
| Thyroid hormone receptor alpha / beta | THRA / THRB (Tox21: TR, no alpha/beta distinction) |
| PPAR gamma | PPARG |

The insulin receptor (IR) is not yet covered.

## Architecture: three independent layers

The layers are never merged, because they measure different things and the assays are not comparable to each other.

| Layer | Folder | Models | Content |
|---|---|---|---|
| ChEMBL IC50 | `models_final/` | 18 | pIC50 regression and classification, by target and assay type (binding `B` or functional `F`) |
| ChEMBL Ki, EC50 and Kd | `models_layers/` | 32 | same targets, each measure in its own model |
| Tox21 qHTS | `models_tox21/` | 11 | classification by assay (environmental chemistry), with viability counter-assays |

### Alert rule

For each molecule-receptor pair, only layers in which the molecule is **inside the applicability
domain** (maximum Tanimoto similarity to the training set greater than or equal to the model's
threshold) and whose model is deployable count. Over these layers:

| Alert | Criterion |
|---|---|
| high (pharmacological + environmental) | 2 or more active layers, with at least one ChEMBL and one Tox21 |
| moderate (majority of layers) | fraction of active layers of at least 0.5, with 3 or more layers assessed |
| low (isolated evidence) | at least 1 active layer |
| no alert | no active layer |

A pair in which no layer is within the domain appears as `-` in the table and as `n.d.` in the plots.

## Installation

Requires Python 3.11.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# check: should print 2180
python -c "import features as F; print(2048 + len(F.DESC_NAMES))"
```

The versions of `rdkit`, `scikit-learn`, and `xgboost` must match those in `requirements.txt`. With a different version
of RDKit, the number of descriptors changes and the predictor fails with `Feature shape mismatch`.
## Usage

```bash
python run_edc.py --input demo_edcs.csv --name-col name --out demo
```

The input CSV needs a `smiles` column (or specify `--smiles-col`). Each run writes to
`results/AAAAMMDD_HHMMSS_<name>/`:

- `<name>_summary.csv`: alert and layer count per molecule and receptor
- `<name>_detail.csv`: result of each model, with domain, probability or predicted value, and conformal interval
- `<name>_rejected.csv`: molecules rejected during standardization, when applicable
- `<name>_dashboard.png` and `<name>_active_fraction.png`: alert plots
- `input.csv` and `README.txt`: copy of the input, command, date, and library versions

It is also possible to use the scripts separately:

```bash
python predict_edc2.py --input demo_edcs.csv --name-col name --out demo   # prediction only
python plot_edc_results.py --out demo                                     # plots only, from demo_*.csv
python predict_edc2.py --smiles "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1"           # a single molecule
```

`panel_ref.csv` is a reference panel with known classes and `demo_edcs.csv` is a short example.

## Repository contents

| Path | Function |
|---|---|
| `run_edc.py` | runs the prediction, generates the plots, and organizes the output in `results/` |
| `predict_edc2.py` | multilayer consensus predictor |
| `plot_edc_results.py` | alert and active-fraction plots |
| `features.py`, `curation.py` | featurization (count-based ECFP4 and RDKit descriptors) and structure standardization |
| `extract_chembl.py`, `extract_bindingdb.py` | extraction of activities from the ChEMBL 37 and BindingDB dumps |
| `pipeline_curate.py`, `pipeline_layers.py`, `pipeline_tox21.py` | curation of the IC50, Ki/Kd/EC50, and Tox21 datasets |
| `modeling.py`, `modeling_layers.py`, `modeling_tox21.py` | training and validation per layer |
| `finalize_models.py`, `conformal.py` | final models and conformal prediction with optimized decision threshold |
| `datasets_layers_by_endpoint/`, `datasets_tox21_by_assay/` | curated Ki/EC50/Kd and Tox21 datasets used in training |
| `models_final/`, `models_layers/`, `models_tox21/` | trained models (about 110 MB) |

Not included: the raw ChEMBL and BindingDB dumps, needed to redo the extraction, and the
curated IC50-layer datasets.
## Validation protocol

- External split by Bemis-Murcko scaffold (80/20) and model selection via 5-fold cross-validation on the training set.
- Final model refit on 100% of the curated data, with metrics coming from the external scaffold-based test.
- Y-scrambling and applicability domain based on Tanimoto similarity.
- Inductive conformal prediction with nominal coverage of 80%, and decision threshold chosen via cross-validation (maximum MCC) instead of 0.5.
- Curation: valid structures, standardized SMILES, removal of salts and solvents, duplicate detection, and stereochemistry checking. Measure types (IC50, Ki, Kd, EC50) and assay types are never mixed.

## Limitations

- Computational screening. An alert indicates plausibility and does not replace experimental assay or regulatory assessment.
- Molecules outside the applicability domain do not enter the alert rule, even when the model returns a high probability. Compounds from classes poorly represented in training, such as phthalates, may end up with no coverage.
- Mechanisms that do not go through the modeled receptors are not detected. For example, effects via aromatase (CYP19A1) or steroidogenesis would require additional endpoints.
- The Tox21 layer is environmental chemistry and uses a different measure type (AC50) than the ChEMBL layers.

## Data sources

ChEMBL 37, BindingDB and Tox21. Derived data follow the terms of each source. ChEMBL is distributed under CC BY-SA 3.0. Check the terms of the other sources before redistributing the datasets.

## License

To be determined.
