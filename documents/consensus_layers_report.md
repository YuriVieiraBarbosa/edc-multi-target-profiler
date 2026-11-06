# Parallel bioactivity layers and consensus prediction

Extension of the endocrine-disruption QSAR pipeline to overcome the three
limitations identified in the first delivery: scope restricted to IC50, absence
of environmental chemistry, and a binary applicability domain.

---

## 1. Decisions made

| Limitation | Decision | Measured result |
|---|---|---|
| IC50 only — genistein and rosiglitazone were left out | Extract Ki, Kd, and EC50 from ChEMBL 37 as **independent layers**, never merged with IC50 | 17 new datasets, 7,733 compound–endpoint pairs |
| Environmental contaminants outside the domain | Incorporate the complete **Tox21 qHTS** (library of ~10,500 substances) via PubChem FTP file | 12 datasets, 7226 unique compounds |
| "Activity" in qHTS confounds antagonism with cell death | Use the **paired viability counter-assays** and autofluorescence controls | 15% to 55% of antagonism actives reclassified |
| Binary applicability domain, with no per-compound uncertainty | **Inductive conformal prediction** with declared and verified nominal coverage | empirical coverage 0.57–0.89 on novel scaffolds |
| Imbalanced classifiers locked at 0.5 | Optimize the threshold via cross-validation and **keep 0.5 when the optimization does not transfer** | only functional PPARγ improved (MCC 0.82 → 0.91) |
| Single-target prediction, with no integration | **Consensus predictor** with per-layer evidence and an alert score | sensitivity 12/15, specificity 3/3 on the reference panel |

![Layer expansion]({{artifact:art_521608b2-0ad8-4598-8310-682fd8724f91}})

---

## 2. Data layers

Each measure kind is its own endpoint. Assay comparability remains the
criterion governing everything: human SINGLE PROTEIN target, `confidence_score ≥ 8`,
`=` relation, unit convertible to nM, median per InChIKey with replicates
discarded when the range exceeds 1 log unit, and separation between binding
and functional assay.
### 2.1 ChEMBL 37 — Ki, Kd and EC50

15,007 raw records, 10,273 after filters.

| receptor   | measure   | assay    |   records |   final compounds |   pAct_median |   discarded_discordant |
|:-----------|:---------|:----------|------------:|-------------------:|---------------:|---------------------------:|
| PPARγ      | EC50     | binding   |        2356 |               1808 |           5.96 |                         35 |
| PPARγ      | EC50     | functional |        1764 |               1342 |           6.52 |                         51 |
| AR         | Ki       | binding   |        1155 |                855 |           7.21 |                         36 |
| ERα        | EC50     | binding   |         822 |                661 |           7.28 |                         35 |
| ERα        | Ki       | binding   |         648 |                469 |           6.99 |                         22 |
| PPARγ      | Ki       | binding   |         566 |                466 |           6.6  |                         12 |
| ERβ        | Ki       | binding   |         405 |                351 |           6.92 |                          4 |
| ERβ        | EC50     | binding   |         411 |                290 |           6.69 |                         17 |
| AR         | EC50     | functional |         311 |                271 |           7.82 |                          3 |
| ERβ        | EC50     | functional |         350 |                254 |           6.98 |                         19 |
| ERα        | EC50     | functional |         293 |                241 |           6.48 |                          6 |
| AR         | EC50     | binding   |         278 |                232 |           8.01 |                          6 |
| THRβ       | EC50     | binding   |         178 |                134 |           6.89 |                          7 |
| THRα       | EC50     | binding   |         149 |                111 |           6.66 |                          4 |
| THRβ       | Ki       | binding   |         118 |                104 |           8.47 |                          2 |
| THRα       | Ki       | binding   |          88 |                 80 |           6.97 |                          1 |
| PPARγ      | Kd       | binding   |          81 |                 65 |           6.02 |                          1 |
### 2.2 Tox21 — qHTS AC50

Twelve nuclear receptor assays from the PubChem FTP file (block
2283001–2284000), with the paired viability counter-assays and the
autofluorescence controls in HEK293. AC50 is given in µM and is reported by
PubChem even for inactive substances — it is only used when the outcome is *Active*.

| receptor   | mode        | cell_line   |   n_tested |   n_active |   n_cytotox_confounds |   valid_compounds |   frac_actives |
|:-----------|:------------|:------------|-------------:|-----------:|---------------------:|--------------------:|--------------:|
| ERα        | agonist    | HEK293T-BLA |        10496 |        594 |                    0 |                6965 |        0.0441 |
| ERα        | antagonist | HEK293T-BLA |        10496 |       1315 |                  198 |                4871 |        0.1482 |
| ERα        | agonist    | VM7-LUC     |        10496 |       1283 |                    0 |                5158 |        0.1483 |
| ERα        | antagonist | VM7-LUC     |        10496 |       1028 |                  316 |                4817 |        0.0924 |
| AR         | agonist    | HEK293T-BLA |        10496 |        625 |                    0 |                4845 |        0.0747 |
| AR         | antagonist | HEK293T-BLA |        10496 |       1915 |                  761 |                4839 |        0.1469 |
| AR         | agonist    | MDA-kb2-LUC |        10496 |        501 |                    0 |                6580 |        0.0433 |
| AR         | antagonist | MDA-kb2-LUC |        10496 |       1196 |                  569 |                4626 |        0.0841 |
| PPARγ      | agonist    | HEK293-BLA  |        10496 |        296 |                    0 |                5777 |        0.0301 |
| PPARγ      | antagonist | HEK293-BLA  |         9315 |        858 |                  198 |                4522 |        0.1    |
| THR*       | agonist    | GH3-LUC     |        10496 |         35 |                    0 |                6729 |        0.0034 |
| THR*       | antagonist | GH3-LUC     |        10496 |       2651 |                 1459 |                3458 |        0.2224 |

Two curation rules specific to this layer:

1. **Cytotoxicity confounding.** An active compound whose viability AC50 is
   less than or equal to 3× the AC50 of the target assay is not distinguishable
   from cell death. These records receive a missing label — they are **not**
   reclassified as inactive, because the information about the receptor
   simply does not exist in that measurement.
2. **Replicates.** Tox21 tests several samples (SIDs) of the same compound,
   and they disagree. Resolution is by majority vote among the conclusive
   replicates, with ties discarded. Applying this rule improved the ROC-AUC
   on the scaffold test from 0.881 to 0.922 for AR agonist (HEK293T) and from
   0.841 to 0.884 for AR agonist (MDA-kb2) — the gain is label noise removal,
   not overfitting.
### 2.3 Layers are not interchangeable

For compounds in common, the correlation between layers ranges from r = −0.24 to r = 0.90,
with median absolute error of 0.4 to 0.9 log units even in the most
correlated pairs. The most discordant pair — AR EC50 functional versus AR IC50
functional — measures agonism versus antagonism in the same assay type, and is
precisely where merging the data would produce the worst damage.

| layer_A      | layer_B      |   n_common |   r_pearson |   rho_spearman |   bias_A_minus_B |   median_absolute_error |
|:-------------|:-------------|----------:|------------:|---------------:|-----------------:|------------------------:|
| PPARG|EC50|F | PPARG|IC50|B |       400 |       0.809 |          0.812 |             0.01 |                    0.41 |
| PPARG|EC50|B | PPARG|IC50|B |       277 |       0.619 |          0.635 |            -0.24 |                    0.68 |
| AR|IC50|F    | AR|Ki|B      |       265 |       0.568 |          0.557 |            -0.21 |                    0.45 |
| AR|EC50|F    | AR|Ki|B      |       205 |       0.538 |          0.555 |            -0.08 |                    0.52 |
| PPARG|EC50|F | PPARG|Ki|B   |       177 |       0.735 |          0.736 |            -0.23 |                    0.46 |
| AR|IC50|B    | AR|Ki|B      |       176 |       0.567 |          0.636 |            -0.05 |                    0.63 |
| ESR1|EC50|B  | ESR1|IC50|B  |       148 |       0.61  |          0.582 |             0.51 |                    0.91 |
| ESR2|EC50|F  | ESR2|IC50|B  |       134 |       0.862 |          0.861 |             0.24 |                    0.42 |
| AR|EC50|B    | AR|IC50|B    |       118 |       0.763 |          0.734 |            -0.03 |                    0.77 |
| ESR1|IC50|B  | ESR1|Ki|B    |       104 |       0.789 |          0.771 |            -0.46 |                    0.61 |
| PPARG|EC50|F | PPARG|IC50|F |        91 |       0.759 |          0.758 |            -0.11 |                    0.46 |
| ESR1|EC50|F  | ESR1|IC50|B  |        86 |       0.685 |          0.718 |             0.46 |                    0.62 |

Tox21 barely overlaps with ChEMBL: between 10 and 26 compounds in common per
pair. This confirms that it covers complementary, not redundant, chemical space.

---

## 3. Models
### 3.1 Ki / Kd / EC50 Layers

Same protocol as the IC50 models: external split by Bemis-Murcko
scaffold, selection via 5-fold CV, y-scrambling, production refit on 100%
of the data.

| receptor   | measure   | assay     |    n | best_model   |   Q2_cv |   y_scrambling_Q2_mean |   coverage_conformal_scaffold | not_deployable   |
|:-----------|:---------|:----------|-----:|:----------------|--------:|------------------------:|-------------------------------:|:------------------|
| AR         | EC50     | binding   |  232 | XGBoost         |   0.647 |                  -0.365 |                          0.783 | False             |
| AR         | EC50     | functional |  271 | XGBoost         |   0.276 |                  -0.427 |                          0.5   | True              |
| AR         | Ki       | binding   |  855 | RandomForest    |   0.649 |                  -0.173 |                          0.684 | False             |
| ERα        | EC50     | binding   |  661 | RandomForest    |   0.66  |                  -0.104 |                          0.636 | False             |
| ERα        | EC50     | functional |  241 | Ridge           |   0.474 |                  -0.58  |                          0.812 | False             |
| ERα        | Ki       | binding   |  469 | XGBoost         |   0.662 |                  -0.348 |                          0.574 | False             |
| ERβ        | EC50     | binding   |  290 | XGBoost         |   0.402 |                  -0.37  |                          0.672 | False             |
| ERβ        | EC50     | functional |  254 | XGBoost         |   0.661 |                  -0.347 |                          0.647 | False             |
| ERβ        | Ki       | binding   |  351 | RandomForest    |   0.688 |                  -0.18  |                          0.771 | False             |
| PPARγ      | EC50     | binding   | 1808 | XGBoost         |   0.646 |                  -0.244 |                          0.682 | False             |
| PPARγ      | EC50     | functional | 1342 | XGBoost         |   0.551 |                  -0.299 |                          0.698 | False             |
| PPARγ      | Kd       | binding   |   65 | Ridge           |   0.727 |                  -0.502 |                          1     | True              |
| PPARγ      | Ki       | binding   |  466 | XGBoost         |   0.609 |                  -0.363 |                          0.656 | False             |
| THRα       | EC50     | binding   |  111 | RandomForest    |   0.471 |                  -0.144 |                          0.636 | True              |
| THRα       | Ki       | binding   |   80 | SVR             |   0.361 |                  -0.338 |                          0.812 | True              |
| THRβ       | EC50     | binding   |  134 | XGBoost         |   0.678 |                  -0.327 |                          0.815 | True              |
| THRβ       | Ki       | binding   |  104 | XGBoost         |   0.54  |                  -0.397 |                          0.571 | True              |
### 3.2 Tox21 Classifiers

| endpoint                       | model_name   |   ROC_AUC_cv |   ROC_AUC_test |   PR_AUC_test |   MCC_test_threshold |   threshold |
|:-------------------------------|:-------------|-------------:|----------------:|---------------:|-------------------:|---------:|
| TOX21_AR_agonist_HEK293T       | RandomForest |        0.844 |           0.922 |          0.773 |              0.692 |     0.4  |
| TOX21_AR_agonist_MDA           | RandomForest |        0.788 |           0.884 |          0.788 |              0.788 |     0.36 |
| TOX21_AR_antagonist_HEK293T    | XGBoost      |        0.834 |           0.782 |          0.274 |              0.274 |     0.24 |
| TOX21_AR_antagonist_MDA        | XGBoost      |        0.84  |           0.832 |          0.375 |              0.343 |     0.15 |
| TOX21_ESR1_agonist_HEK293T     | RandomForest |        0.882 |           0.759 |          0.217 |              0.117 |     0.26 |
| TOX21_ESR1_agonist_VM7         | RandomForest |        0.732 |           0.814 |          0.661 |              0.452 |     0.35 |
| TOX21_ESR1_antagonist_HEK293T  | XGBoost      |        0.856 |           0.805 |          0.426 |              0.361 |     0.11 |
| TOX21_ESR1_antagonist_VM7      | RandomForest |        0.862 |           0.833 |          0.444 |              0.404 |     0.2  |
| TOX21_PPARG_agonist_HEK293     | RandomForest |        0.834 |           0.796 |          0.22  |              0     |     0.4  |
| TOX21_PPARG_antagonist_HEK293  | XGBoost      |        0.836 |           0.829 |          0.465 |              0.366 |     0.13 |
| TOX21_TR_antagonist_GH3        | XGBoost      |        0.8   |           0.72  |          0.235 |              0.253 |     0.28 |

The TR agonist (GH3) assay was discarded: 35 actives among 10,496 substances do not
support a model.
### 3.3 Conformal calibration and operating point

| endpoint   | model_name   |   q_absolute |   coverage_scaffold | conformal_mode   |   production_threshold |   threshold_gain |
|:-----------|:-------------|-------------:|---------------------:|:-----------------|------------------:|---------------:|
| AR_B       | XGBoost      |        0.801 |                0.673 | fixed            |              0.5  |         -0.041 |
| AR_F       | RandomForest |        0.764 |                0.79  | fixed            |              0.5  |          0     |
| ESR1_B     | XGBoost      |        0.837 |                0.772 | fixed            |              0.5  |         -0.063 |
| ESR1_F     | RandomForest |        1.234 |                0.645 | fixed            |              0.5  |          0     |
| ESR2_B     | XGBoost      |        0.925 |                0.565 | fixed            |              0.5  |         -0.06  |
| PPARG_B    | XGBoost      |        0.82  |                0.718 | fixed            |              0.5  |         -0.035 |
| PPARG_F    | Ridge        |        1.025 |                0.739 | fixed            |              0.57 |          0.092 |
| THRA_B     | RandomForest |        0.592 |                0.705 | normalized       |              0.5  |         -0.003 |
| THRB_B     | XGBoost      |        0.715 |                0.886 | fixed            |              0.5  |          0     |

**Conformal coverage falls below the nominal level.** With α = 0.20, the
empirical coverage on the scaffold-based test set ranges from 0.57 to 0.89 (median 0.72), not 0.80. This
is not an implementation defect: the conformal guarantee holds under exchangeability,
and a scaffold-based split breaks exchangeability on purpose — that is the reason
for using it. The honest number for interpreting a prediction on new chemistry is the
coverage measured on new scaffolds, recorded in
`coverage_scaffold` within each model, not the nominal α.

Threshold optimization via CV **did not transfer**: of the 9 IC50 endpoints, the MCC
on the scaffold-based test worsened in 5 (AR binding, ERα binding, ERβ binding, PPARγ
binding, THRα binding), remained unchanged in 3 (AR functional, ERα functional, THRβ
binding), and improved in 1. The optimized threshold was adopted only for PPARγ functional. Three endpoints are flagged as non-deployable
in the model file itself: AR functional, ERα functional, and THRα binding.

---
## 4. Consensus prediction

`predict_edc2.py` evaluates each molecule across all available layers and returns
two tables: a detailed one, with one row per layer (predicted value, conformal
interval, probability, similarity to the nearest training neighbor, and domain
flag), and a summary per receptor with the alert score.

| Level | Rule |
|---|---|
| high | at least 2 active layers, with at least one pharmacological (ChEMBL) and one environmental (Tox21) |
| moderate | at least 3 layers assessed and majority active |
| low | isolated evidence, 1 active layer |
| no alert | no active layer within the domain |

Only layers in which the molecule is within the applicability domain **and**
whose model is not flagged as non-deployable are included in the score.

### Reference panel

19 input compounds: 9 recognized endocrine disruptors, 6 mechanism references
(ERα agonist, PPARγ agonist, AR antagonist, phytoestrogen,
estrogenic mycotoxin, and antimicrobial), 3 negative controls, and 1 non-EDC drug.

**18 were assessed.** Glycine is rejected during standardization due to molar mass
below 100 Da — the pipeline's MW filter, applied to all sources. The
predictor now writes the rejected molecules to a `_rejected.csv` file and
reports them in the output, rather than simply omitting them.

**Sensitivity 12/15, specificity 3/3** (among the 18 assessed; of the 4
non-EDC controls in the panel, glycine was not evaluated). The two gaps that motivated this
phase were closed: genistein now reaches high ERα (via the Ki and EC50 layers) and
rosiglitazone reaches high PPARγ. DDT, methoxychlor, vinclozolin, and triclosan, previously
outside the domain, are now recovered by the Tox21 layer.

The three false negatives — atrazine, DBP, and DEHP — deserve careful reading:
all three are **measured as inactive** in the Tox21 assays themselves for these five
receptors. Atrazine acts via aromatase induction, DEHP via the MEHP metabolite, and
DBP via interference with fetal steroidogenesis. These are not model failures; they are
mechanisms outside the modeled receptor's scope. Expanding coverage for them
requires aromatase (CYP19A1) and steroidogenesis (H295R) endpoints, not
better models of the same receptors.

---
## 5. Remaining limitations

1. **ERβ has no environmental layer.** Tox21 does not include an ERβ assay. The
   α/β selectivity remains supported only by pharmacological data.
2. **THR has no subtype resolution.** The Tox21 TR assay is in rat GH3 cells
   and does not distinguish α from β; in the tables it appears as `THR*` and is never
   added to the human endpoints THRalpha and THRbeta.
3. **PPARgamma agonist in Tox21 has a prevalence of 3%.** The MCC on the scaffold
   test is 0, despite the ROC-AUC of 0.80: the model_name ranks, but does not classify.
   Use the probability as a ranking, not as a decision.
4. **Conformal coverage is lower than nominal for novel chemistry** (section 3.3).
5. **Cortellis remains outside.** Clarivate's licensed product, with no
   public API. A CSV export enters the same provenance structure.
6. **Nothing here replaces an in vitro assay** nor constitutes a regulatory risk
   assessment. The alert score is a prioritization tool.

---

## 6. Reproducibility

| File | Content |
|---|---|
| `pipeline_layers.py` | extraction and curation of the Ki / Kd / EC50 layers from ChEMBL 37 |
| `pipeline_tox21.py` | Tox21 curation: cytotoxicity, autofluorescence, replicates |
| `modeling_layers.py` | models for the Ki / Kd / EC50 layers |
| `modeling_tox21.py` | Tox21 classifiers with optimized threshold |
| `conformal.py` | conformal calibration and operating point of the IC50 models |
| `predict_edc2.py` | multi-layer consensus predictor |

Usage:

```bash
python predict_edc2.py --input molecules.csv --smiles-col smiles --name-col name
```

Environment: Python 3.11, RDKit 2026.03.6, scikit-learn 1.9.1, XGBoost. Fixed
seed (42) across all splits and models.
