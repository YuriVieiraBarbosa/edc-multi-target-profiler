# Engineering the EDC-QSAR pipeline: data, models, and consensus

This document describes what the code in `~/Downloads/pipeline_edc_qsar` does, in the order in which the data flow. The statements come from reading the scripts (`curation.py`, `features.py`, `pipeline_curate.py`, `pipeline_layers.py`, `pipeline_tox21.py`, `modeling.py`, `modeling_layers.py`, `modeling_tox21.py`, `finalize_models.py`, `conformal.py`, `predict_edc2.py`), the `README_TEXT_*.csv` files that accompany the models, the 61 `.joblib` files loaded, and two runs in `results/`. The curated datasets for the IC50 layer are not in the folder; for IC50 the tables use the `README_TEXT` files and the metadata recorded in the models.

## 1. Overview

```
ChEMBL 37 (SQLite) --+-- IC50 ---------> pipeline_curate.py -> dataset_<target>_<B|F>.csv --> modeling.py -> finalize_models.py -> conformal.py --> models_final/  (18)
                     +-- Ki, Kd, EC50 -> pipeline_layers.py -> dataset_<target>_<measure>_<B|F>.csv -> modeling_layers.py --------------------------> models_layers/ (32)
Tox21 (PubChem qHTS) -> pipeline_tox21.py -> dataset_TOX21_<receptor>_<mode>_<cell_line>.csv -> modeling_tox21.py ----------------------------------> models_tox21/ (11)
BindingDB ----------> pipeline_curate.py -> external_bindingdb.csv (external validation of the IC50 models)

61 models --> predict_edc2.py --> 1 binary vote per model --> alert per molecule x receptor
```

Architecture rule: each measure type (IC50, Ki, Kd, EC50, AC50) and each assay type (binding `B`, functional `F`) has its own dataset and its own pair of models. Nothing is merged before the consensus.

## 2. Data
### 2.1 Why the measures should not be mixed

| Measure | What it measures | Assay dependency |
|---|---|---|
| IC50 | concentration that inhibits 50% of a response | high: varies with ligand/substrate concentration and with assay format (relates to Ki via Cheng-Prusoff in competitive binding assays) |
| Ki / Kd | equilibrium dissociation constant | low in principle (thermodynamic quantity) |
| EC50 | concentration of half-maximal effect in a functional assay | high: depends on the cellular system, the reporter, and signal amplification |
| AC50 (Tox21) | concentration of half-maximal activity, derived from a Hill fit on a fixed concentration grid | covers activation and inhibition; it is also reported for inactives; in antagonist mode it suffers confounding by cytotoxicity |

The numerical scales coincide: pIC50 = 9 − log10(nM) and pAC50 = 6 − log10(µM) are both −log10(M). This is why a threshold of "6" seems comparable between layers, but the biological meaning is not the same. This is the reason for keeping the layers separate and only integrating them at the consensus stage.
### 2.2 ChEMBL Curation (IC50 in `pipeline_curate.py`; Ki, Kd and EC50 in `pipeline_layers.py`)

The two scripts apply the same funnel. The only difference is the `standard_type` queried.

1. Targets: AR (P10275), ERα (P03372), ERβ (Q92731), PPARγ (P37231), THRα (P10827) and THRβ (P10828), restricted to `SINGLE PROTEIN` from *Homo sapiens*.
2. `confidence_score >= 8` (direct protein target).
3. `standard_relation == '='` (discards `<`, `>`, `~`).
4. Unit convertible to nM (nM, µM, mM, M, pM).
5. Value between 0.001 nM and 10⁷ nM (pActivity from 2 to 12).
6. Discards records with a problematic `data_validity_comment` (outside the typical range, non-standard unit, possible transcription or author error, potential duplicate).
7. Requires a defined chemical structure (`structure_type == 'MOL'`).
8. Structural standardization (`curation.standardize`): cleanup, parent fragment (removes salts and solvents), normalization, charge neutralization, canonical tautomer, InChIKey. Rejects invalid SMILES, non-organic elements, fewer than 6 heavy atoms, and molar mass outside 100–1000 Da, with the reason recorded.
9. Conversion: `pAct = 9 − log10(value in nM)`. The column is called `pIC50` in the IC50 dataset and `pAct` in the Ki/Kd/EC50 layers.
10. Replicates: groups by InChIKey, uses the **median**, and discards the compound if the spread between replicates exceeds 1.0 log unit.
11. Splits by (target, measure type, `assay_type`). Only `B` and `F` are included, with a minimum of 50 records per group.

Each final row retains the provenance: ChEMBL ID, name, original and standardized SMILES, type/value/unit, pChEMBL, assay, target, organism, document, DOI, PMID, and the lists of all assays and documents supporting the compound.

The `assay_mode` (agonist, antagonist, undefined) is annotated by regular expression on the assay description and remains only as a column. No model uses it. For this reason, the ChEMBL models do not distinguish agonism from antagonism.

BindingDB is used only as external validation of the IC50 models: the compounds kept are those whose InChIKey does not exist in ChEMBL for the same target (817 compounds according to the project record). The script that computes the metrics for this validation is not in the folder.
### 2.3 Tox21 / AC50 (`pipeline_tox21.py`)

The Tox21 layer is not an AC50 regression. What enters training is the active/inactive label; the pAC50 is calculated and stored, but not used.

1. Source: 12 PubChem qHTS assays fixed in the `ASSAYS` table (AR and ESR1 in two cell lines each, PPARγ and TR in one, all in agonist and antagonist mode). Each assay is a dataset.
2. Label: `Active` = 1, `Inactive` = 0, `Inconclusive` = empty (excluded from training).
3. Censored points: AC50 exists for all substances, including inactive ones. `pAC50 = 6 − log10(AC50 in µM)` is only filled in when the outcome is `Active`.
4. Cytotoxicity (antagonist modes): if the viability counter-assay is active and its AC50 is less than or equal to 3 times the AC50 of the target assay, the "active" outcome is not distinguishable from cell death. It becomes empty and is excluded from training (it does not become inactive). The 3-fold ratio is a pipeline heuristic, without formal calibration. Agonist assays have no viability counter-assay.
5. Autofluorescence: in BLA assays, an "active" whose compound is active in any of the three HEK293 autofluorescence assays (blue, green, red) is also excluded.
6. Replicates: Tox21 tests multiple SIDs per compound. It is resolved by majority vote among the conclusive replicates, and ties are discarded.
7. Structure: the same standardization as ChEMBL; the compound key is the InChIKey.
8. Training skips assays with fewer than 30 examples in the minority class. This is why `TOX21_TR_agonist_GH3` (27 actives) has no model, and TR only exists in antagonist mode.
### 2.4 Points of attention in the data layer

1. **Stereochemistry** is recorded (`undefined_stereo` counted undefined chiral centers), but it is not used as a filter. Enantiomers have different InChIKeys and are entered as distinct compounds.
2. **Ki and EC50 share compounds.** For the same receptor, 205 of the 271 compounds in `AR_EC50_F` are also in `AR_Ki_B`; for PPARgamma, 177 of the 1,342 compounds in `PPARG_EC50_F` are in `PPARG_Ki_B`. The measurements are distinct, but the structures repeat, so the votes from these layers are not independent.
3. **ChEMBL and Tox21 barely overlap.** Only 0.4% to 4% of the compounds in each ChEMBL dataset (Ki/EC50/Kd) appear in the Tox21 assays for the same receptor (maximum 13.8%, `PPARG_Kd_B`). The two chemical worlds are different: medicinal chemistry ligands versus environmental chemicals.
4. **Discrepancy to investigate.** The CSVs in `datasets_tox21_by_assay/` do not exactly match `README_models_tox21.csv`. For `TOX21_AR_agonist_HEK293T` the CSV has 4,890 compounds and 424 actives (8.7%), while the README records n = 4,844 and an active fraction of 0.075 (about 363 actives). I did not verify the cause; it is consistent with an earlier version of the curation.

## 3. Models

### 3.1 Features (`features.py`)

2,180 columns: ECFP4 by count (radius 2, 2,048 positions) plus 132 2D descriptors from RDKit (all except the `fr_*` ones). Descriptors with NaN or infinity become 0. The binary ECFP4 of the same molecule is used only for the applicability domain. The number 2,180 depends on the RDKit version: with a different version the predictor fails with `Feature shape mismatch`.
### 3.2 Training protocol (same across the three modeling scripts)

1. External 80/20 split by Bemis-Murcko scaffold (`features.scaffold_split`).
2. On the training set, for each algorithm in the "zoo", a small hyperparameter grid with 5-fold CV.
3. Selection of the best algorithm: CV Q² for regression, CV ROC-AUC for classification (Tox21 also uses ROC-AUC).
4. Evaluation on the scaffold-based external test (R², RMSE, and MAE; ROC-AUC, PR-AUC, balanced accuracy, and MCC).
5. y-scrambling on the best regressor (Q² with shuffled labels should be close to zero or negative).
6. Applicability domain and conformal calibration (3.4).
7. Refitting of the best model on 100% of the curated data. This is the production model. Metrics still come from the protocol with the external test.

| Task | Algorithms tested |
|---|---|
| Regression | Ridge, SVR (RBF), RandomForest (500 trees), XGBoost (600 trees) |
| ChEMBL Classification | Logistic, SVC (RBF), RandomForest (balanced), XGBoost |
| Tox21 Classification | Logistic (balanced), RandomForest (balanced), XGBoost |

### 3.3 Regression and classification of the same endpoint

- Regression predicts `pAct` (continuous).
- Classification predicts `pAct >= 6.0` (IC50, Ki, Kd, or EC50 up to 1 µM). The label comes from the same measured value, so **classification does not use new data**: it is the same set of compounds, with the label binarized.
- A classification is only trained if the fraction of actives is between 5% and 95% and the minority class has at least 10 compounds in training. `AR_EC50_F` (96% active) and `THRB_Ki_B` (96% active) have no classifier; this is the reason there are 32 models in `models_layers/` and not 34.
- In the Tox21 layer only classification exists, and the label is the depositor's active/inactive call, not a potency cutoff.
### 3.4 Domain, conformal, and threshold

- **Applicability domain:** the threshold is the 5th percentile of Tanimoto similarity to the nearest neighbor within the training set itself. A new molecule is within the domain if its maximum similarity to the training set is greater than or equal to this threshold. Values: 0.29 to 0.63 for the ChEMBL models and 0.265 to 0.292 for the Tox21 ones. By construction, 5% of the training set itself falls outside, and 0.27 is a low similarity, so the Tox21 domain is permissive.
- **Inductive conformal (regression), α = 0.20:** splits the training set into fitting (75%) and calibration (25%), takes the quantile of the absolute residuals, and records `interval_lower` and `interval_upper`. The nominal coverage is 80%; the coverage measured on new scaffolds is recorded in `coverage_scaffold`. The quantile is calibrated on a model fitted with only 75% of the training set, while the production model is refitted with 100%.
- **Classification decision threshold:** in `models_final`, `conformal.py` computes the threshold that maximizes the MCC in CV, and it was only adopted where it improved the external test (only `PPARG_F`, threshold 0.57; for the rest it stays at 0.5). In `models_layers` the threshold is fixed at 0.5. For the Tox21 models the optimized threshold is always adopted, and it ranges from 0.11 to 0.40 because the classes are highly imbalanced.
- **`not_deployable`:** regression in `models_layers` if n < 200 or CV Q² < 0.3; classification in `models_layers` if n < 200; in `models_final` the flag comes from `README_calibration.csv`. These models are recorded, but the consensus ignores them.
### 3.5 How the scaffold split behaves (verified)

`scaffold_split` iterates over the scaffold groups from largest to smallest and sends each group to the test set if it fits within the 20% budget. The docstring says the opposite ("larger scaffolds go to training"). I ran the function on the 17 Ki/Kd/EC50 datasets in the folder:

- The test set ends up with the most populous chemical families and the training set with the rest. Examples: `AR_Ki_B` has 171 test compounds coming from only 4 scaffolds; `PPARG_EC50_B` has 362 compounds coming from 13 scaffolds. In all datasets the test set has between 2 and 15 scaffolds.
- The average group size per scaffold ranges from 2.7 to 42.8 compounds in the test set (median 16.5) and from 1.5 to 4.3 in the training set. Among the 5 largest scaffolds in each dataset, 1 to 5 are in the test set.

Consequence: the external test set is an extrapolation test for the dominant chemical families, dependent on few scaffolds and pessimistic. This is consistent with the drop in R² between CV and test (CV is random within the training set) and with conformal coverage below nominal. Test metrics should be read as a lower bound, with high variance. The production model, refit on 100%, contains these families. I did not change anything; it remains to decide what the intention was and, if necessary, swap the rule (large groups to training) and reevaluate. I only checked the Ki/Kd/EC50 datasets; the IC50 ones use the same function, but their CSVs are not in the folder.
### 3.6 Models and performance recorded in the README_TEXT

**ChEMBL IC50 (`models_final/`, 9 endpoints, 18 models).** Conformal coverage on novel scaffolds; MCC and balanced accuracy on the scaffold test with threshold 0.5.

| endpoint | n | reg. | clf. | conformal coverage (scaffold) | MCC test | BalAcc test | threshold | not deployable |
|---|---|---|---|---|---|---|---|---|
| AR_B | 1807 | XGBoost | XGBoost | 0.67 | 0.39 | 0.69 | 0.50 | no |
| AR_F | 618 | RandomForest | RandomForest | 0.79 | 0.00 | 0.50 | 0.50 | yes |
| ESR1_B | 2700 | XGBoost | XGBoost | 0.77 | 0.48 | 0.75 | 0.50 | no |
| ESR1_F | 156 | RandomForest | RandomForest | 0.65 | 0.83 | 0.88 | 0.50 | yes |
| ESR2_B | 1548 | XGBoost | RandomForest | 0.56 | 0.39 | 0.75 | 0.50 | no |
| PPARG_B | 1437 | XGBoost | Logistic | 0.72 | 0.71 | 0.81 | 0.50 | no |
| PPARG_F | 228 | Ridge | SVC | 0.74 | 0.82 | 0.91 | 0.57 | no |
| THRA_B | 220 | RandomForest | RandomForest | 0.70 | 0.28 | 0.58 | 0.50 | yes |
| THRB_B | 396 | XGBoost | XGBoost | 0.89 | 0.49 | 0.78 | 0.50 | no |

**ChEMBL Ki, Kd and EC50 (`models_layers/`, 17 endpoints, 32 models).** Q² CV of the best regressor and mean Q² with scrambled labels.

| endpoint | n | n_scaffolds | active frac. | regressor | Q² CV | Q² scrambling | cov. conformal | threshold AD | classifier | not deployable |
|---|---|---|---|---|---|---|---|---|---|---|
| AR_EC50_B | 232 | 78 | 0.77 | XGBoost | 0.65 | -0.36 | 0.78 | 0.40 | XGBoost | no |
| AR_EC50_F | 271 | 52 | 0.96 | XGBoost | 0.28 | -0.43 | 0.50 | 0.53 | no classifier | yes |
| AR_Ki_B | 855 | 238 | 0.84 | RandomForest | 0.65 | -0.17 | 0.68 | 0.58 | RandomForest | no |
| ESR1_EC50_B | 661 | 277 | 0.75 | RandomForest | 0.66 | -0.10 | 0.64 | 0.29 | RandomForest | no |
| ESR1_EC50_F | 241 | 71 | 0.73 | Ridge | 0.47 | -0.58 | 0.81 | 0.42 | Logistic | no |
| ESR1_Ki_B | 469 | 188 | 0.82 | XGBoost | 0.66 | -0.35 | 0.57 | 0.48 | SVC | no |
| ESR2_EC50_B | 290 | 98 | 0.74 | XGBoost | 0.40 | -0.37 | 0.67 | 0.41 | RandomForest | no |
| ESR2_EC50_F | 254 | 76 | 0.75 | XGBoost | 0.66 | -0.35 | 0.65 | 0.44 | RandomForest | no |
| ESR2_Ki_B | 351 | 135 | 0.76 | RandomForest | 0.69 | -0.18 | 0.77 | 0.49 | RandomForest | no |
| PPARG_EC50_B | 1808 | 672 | 0.49 | XGBoost | 0.65 | -0.24 | 0.68 | 0.58 | XGBoost | no |
| PPARG_EC50_F | 1342 | 576 | 0.70 | XGBoost | 0.55 | -0.30 | 0.70 | 0.63 | RandomForest | no |
| PPARG_Kd_B | 65 | 36 | 0.52 | Ridge | 0.73 | -0.50 | 1.00 | 0.30 | Logistic | yes |
| PPARG_Ki_B | 466 | 218 | 0.68 | XGBoost | 0.61 | -0.36 | 0.66 | 0.50 | RandomForest | no |
| THRA_EC50_B | 111 | 41 | 0.72 | RandomForest | 0.47 | -0.14 | 0.64 | 0.62 | SVC | yes |
| THRA_Ki_B | 80 | 28 | 0.82 | SVR | 0.36 | -0.34 | 0.81 | 0.54 | XGBoost | yes |
| THRB_EC50_B | 134 | 46 | 0.86 | XGBoost | 0.68 | -0.33 | 0.81 | 0.62 | RandomForest | yes |
| THRB_Ki_B | 104 | 44 | 0.96 | XGBoost | 0.54 | -0.40 | 0.57 | 0.57 | no classifier | yes |
**Tox21 (`models_tox21/`, 11 classification models).**

| endpoint | n | active frac. | n_scaffolds | model_name | threshold | AD threshold |
|---|---|---|---|---|---|---|
| TOX21_AR_agonist_HEK293T | 4844 | 0.075 | 1462 | RandomForest | 0.40 | 0.278 |
| TOX21_AR_agonist_MDA | 6579 | 0.043 | 2008 | RandomForest | 0.36 | 0.289 |
| TOX21_AR_antagonist_HEK293T | 4838 | 0.147 | 1600 | XGBoost | 0.24 | 0.288 |
| TOX21_AR_antagonist_MDA | 4625 | 0.084 | 1502 | XGBoost | 0.15 | 0.274 |
| TOX21_ESR1_agonist_HEK293T | 6964 | 0.044 | 2182 | RandomForest | 0.26 | 0.292 |
| TOX21_ESR1_agonist_VM7 | 5157 | 0.148 | 1524 | RandomForest | 0.35 | 0.286 |
| TOX21_ESR1_antagonist_HEK293T | 4870 | 0.148 | 1642 | XGBoost | 0.11 | 0.273 |
| TOX21_ESR1_antagonist_VM7 | 4816 | 0.092 | 1601 | RandomForest | 0.20 | 0.273 |
| TOX21_PPARG_agonist_HEK293 | 5776 | 0.030 | 1750 | RandomForest | 0.40 | 0.286 |
| TOX21_PPARG_antagonist_HEK293 | 4522 | 0.100 | 1468 | XGBoost | 0.13 | 0.277 |
| TOX21_TR_antagonist_GH3 | 3457 | 0.222 | 1237 | XGBoost | 0.28 | 0.265 |
## 4. Consensus (`predict_edc2.py`)

### 4.1 From model to vote

Each of the 61 models becomes one row per molecule in `_detail.csv`, with a binary vote `active`:

| Model | Vote "active" |
|---|---|
| Regression | `predicted_value >= 6.0` (the conformal interval is recorded, but does not enter the vote) |
| Classification | `prob_active >= threshold` of the model |

Before counting, each row receives `in_domain` (maximum similarity to the model's own training set greater than or equal to the domain threshold) and `not_deployable`. Only rows with `counted = in_domain & ~not_deployable` count.

### 4.2 From vote to alert

For each molecule-receptor pair, over the rows that count:

- `layers_assessed` = number of distinct models (regression and classification of the same endpoint count separately);
- `layers_active` = active votes; `pharmacological_evidence` = active votes from ChEMBL models; `tox21_evidence` = active votes from Tox21 models.

| Alert | Criterion |
|---|---|
| high | `layers_active >= 2`, with `pharmacological_evidence >= 1` and `tox21_evidence >= 1` |
| moderate | `active_fraction >= 0.5` and `layers_assessed >= 3` |
| low | `layers_active >= 1` |
| no alert | no active layer |

The rules are evaluated in this order, and the first one that applies holds.

### 4.3 Structural limits of the alert by receptor

The number of possible layers depends on how many deployable models exist per receptor. The Tox21 models for TR fall under the label `TR`, separate from THRA and THRB, and there is no Tox21 assay for ERβ. Calculated from the 61 loaded models:

| receptor (in the summary) | deployable ChEMBL models | Tox21 models | possible layers | maximum attainable alert |
|---|---|---|---|---|
| AR | 6 | 4 | 10 | high |
| ESR1 | 8 | 4 | 12 | high |
| ESR2 | 8 | 0 | 8 | moderate |
| PPARG | 10 | 2 | 12 | high |
| THRA | 0 | 0 | 0 | none (does not appear in the summary) |
| THRB | 2 | 0 | 2 | low |
| TR | 0 | 1 | 1 | low |

Consequences: "high" is only attainable for AR, ESR1, and PPARG. ESR2 reaches a maximum of "moderate". THRβ reaches a maximum of "low" (two deployable models, both from the same IC50 endpoint). TR, with a single model, reaches a maximum of "low". THRα has no deployable ChEMBL model and does not appear in the summary. In the two runs in `results/` the observed levels respect these ceilings. Therefore "no alert" or "low" in ERβ, THR, and TR does not carry the same weight as in AR, ESR1, and PPARG.
### 4.4 Consensus points of attention

1. **Double counting of the same endpoint.** Regression and classification of an endpoint use the same compounds and the same features, but count as two layers. In the `demo` and `test` runs there are 80 and 185 layers assessed against 67 and 147 distinct endpoints. Impact measured in these two runs: collapsing to one vote per endpoint (active only if both agree) changes 1 of 22 and 0 of 46 alerts, so the effect is small, but the mechanism exists (for example, Tamoxifen in ESR2 goes from "moderate" to "low").
2. **"Active" has different definitions per layer.** In the ChEMBL layers it means potency up to 1 µM. In Tox21 it means an active curve in the tested range, with median pAC50 between 4.4 and 5.9 across the 11 assays modeled. In the folder's CSVs (see the divergence in 2.4), only about 15% of Tox21 actives would have pAC50 greater than or equal to 6 (from 7% to 48% per assay). Furthermore, the Tox21 probability thresholds (0.11 to 0.40) are much more permissive than ChEMBL's 0.5. A Tox21 active vote and a ChEMBL active vote are not equivalent.
3. **No direction of effect.** The alert sums Tox21 agonist and antagonist votes and ChEMBL votes without mode. It means "there is evidence of interaction with the receptor", not agonism or antagonism.
4. **Equally weighted votes, with different performances.** Outside the training set, classification generalizes better than regression (ROC-AUC from 0.76 to 0.94 in 8 of 9 IC50 endpoints, versus test R² from −0.26 to 0.66). The regression vote weighs the same as the classification vote.
5. **Permissive domain in Tox21** (minimum similarity of about 0.27) and `coverage_scaffold` below nominal (reading in 3.5). This reduces confidence in votes for molecules distant from the training set.

### 4.5 How to read a result

1. See `alert` in `_summary.csv` along with `layers_assessed` and the receptor cap (section 4.3).
2. Open `_detail.csv`, filter `counted = True` and check which models voted. Use `detail` to see the Tox21 mode (agonist or antagonist) and `tanimoto_nn` to see the distance to the training set.
3. Treat `predicted_value` and the conformal interval as ranking and order of magnitude, and `prob_active` as the main indicator of activity.
4. Absence of an alert is not evidence of absence of effect: mechanisms outside the modeled receptors (aromatase, steroidogenesis) and receptors without coverage (IR) are not detected.
## 5. Open decisions

1. Double counting: keep it, use one vote per endpoint, or use only the classification as the vote.
2. `scaffold_split` rule: keep the current one or adopt large groups in training, and reassess the metrics.
3. Direction of effect: document it in the README, add a column with active Tox21 modes to `_summary.csv`, or model by assay and mode (section 9 of the report).
4. Discrepancy between the Tox21 CSVs and the README (item 4 of section 2.4).
