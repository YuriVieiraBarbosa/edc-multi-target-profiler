# QSAR Models for Endocrine Disruption — AR, ERα, ERβ, THRα, THRβ and PPARγ

Curated dataset of IC50 values and machine learning models for screening new
compounds with potential interaction with nuclear receptors associated with
endocrine disruption.

---

## 1. Data sources

| Source | Status | IC50 records across the 6 targets |
|---|---|---|
| **ChEMBL 37** (official SQLite dump, SHA256 verified) | primary source | 19,472 |
| **BindingDB 2026-09** (complete TSV dump) | secondary source / external validation | 17,544 |
| **PubChem BioAssay** (SDQ, by UniProt accession) | consulted, no independent contribution | 20,882, 100% deposited by ChEMBL or BindingDB |
| **Cortellis** (Clarivate) | **not accessible** | — |

**About Cortellis.** It is a licensed commercial product, with no public API and
no credential available in this environment. It was left out of the collection.
If institutional access is available, a CSV export can be integrated into the
same provenance structure without changing the pipeline.

**About the ChEMBL REST API.** During execution, the list endpoints returned
HTTP 500 for about 70% of calls. Extraction was migrated to the official
SQLite dump (`chembl_37_sqlite.tar.gz`, 5.76 GB), queried via SQL — a more
complete and reproducible solution than API pagination.

**Targets.** Resolved by UniProt accession, restricted to human `SINGLE PROTEIN`:

| Receptor | ChEMBL ID | UniProt |
|---|---|---|
| Androgen receptor (AR) | CHEMBL1871 | P10275 |
| Estrogen receptor α (ERα/ESR1) | CHEMBL206 | P03372 |
| Estrogen receptor β (ERβ/ESR2) | CHEMBL242 | Q92731 |
| PPAR-γ | CHEMBL235 | P37231 |
| Thyroid hormone receptor α (THRα) | CHEMBL1860 | P10827 |
| Thyroid hormone receptor β (THRβ) | CHEMBL1947 | P10828 |

Protein complexes, families and selectivity groups (e.g., "Estrogen receptor"
as a family, "PPAR alpha/gamma") were **excluded** — the molecular target must
be unambiguous for the endpoint to be interpretable.

---
## 2. Assay comparability

The central requirement was not to mix biologically incomparable measurements.
Filters applied, in order, to records from ChEMBL:

| filter                                                  |   records remaining |
|:--------------------------------------------------------|----------------------:|
| 0. raw IC50 values (human target SINGLE PROTEIN)        |                 19472 |
| 1. Homo sapiens target                                  |                 19472 |
| 2. confidence_score >= 8 (direct protein target)        |                 19472 |
| 3. standard_relation '='                                |                 14575 |
| 4. unit convertible to nM                                |                 14572 |
| 5. value within plausible range (pIC50 2-12)            |                 14568 |
| 6. no problematic data_validity_comment                  |                 14513 |
| 7. with defined chemical structure                      |                 14502 |
| 8. valid standardized structure (salt/solvent removed)  |                 14445 |

In addition:

- **Separation by assay type.** **Binding** assays (`assay_type = B`,
  typically radioligand displacement) and **functional** assays (`assay_type = F`,
  typically reporter gene in cell) generate **independent** datasets. A competitive
  displacement IC50 and a transcriptional inhibition IC50 measure distinct
  phenomena and were never merged.
- **confidence_score ≥ 8** — direct single-protein assignment. All
  19,472 records already met this criterion.
- **Only the `=` relation.** Censored records (`>`, `<`) were discarded:
  4,897 records, 25% of the raw total.
- **Annotated mode of action** (agonist/antagonist/undefined) derived from the
  assay description, available as a column for additional filtering. The
  functional datasets are dominated by one of the modes (functional AR: 98.7%
  antagonist; functional PPARγ: 7.5% antagonist, i.e., predominantly
  agonist) — mixing these two receptors in the same functional model would be
  incorrect, and they remain separate.
- **Replicates.** Median per InChIKey within each endpoint; compounds with an
  inter-replicate range greater than **1 log unit** were discarded
  (747 compounds in total, 518 of them in ERα binding — 16% of that set).
### Structural standardization (RDKit)

`Cleanup` → selection of the largest fragment (removal of salts and solvents) →
normalization of functional groups → charge neutralization → tautomeric
canonicalization → InChIKey as identity identifier. Rejected: invalid
structures, inorganics, elements outside the organic set, fewer than 6 heavy
atoms, and molecular mass outside 100–1000 Da (30 structures in ChEMBL, 66 in
BindingDB). Undefined stereocenters were **flagged** (column
`undefined_stereo`), not removed: 2,153 of the 9,110 compound–target pairs.

---
## 3. Resulting datasets

**7245 unique compounds, 9,110 compound–endpoint pairs, 9 endpoints.**

| receptor   | assay    |   IC50 records |   compounds |   discarded (discordant replicates) |   final compounds |   median pIC50 |   fraction antagonist mode |
|:-----------|:----------|-----------------:|------------:|--------------------------------------:|-------------------:|----------------:|--------------------------:|
| ERα        | binding   |             5637 |        3218 |                                   518 |               2700 |            7    |                     0.413 |
| AR         | binding   |             2687 |        1883 |                                    76 |               1807 |            6.46 |                     0.603 |
| ERβ        | binding   |             2146 |        1608 |                                    60 |               1548 |            6.89 |                     0.501 |
| PPARγ      | binding   |             1665 |        1464 |                                    27 |               1437 |            6.48 |                     0.153 |
| AR         | functional |              888 |         656 |                                    38 |                618 |            6.7  |                     0.987 |
| THRβ       | binding   |              519 |         402 |                                     6 |                396 |            6.14 |                     0.831 |
| PPARγ      | functional |              300 |         231 |                                     3 |                228 |            6.28 |                     0.075 |
| THRα       | binding   |              272 |         221 |                                     1 |                220 |            5.8  |                     0.714 |
| ERα        | functional |              220 |         174 |                                    18 |                156 |            6.38 |                     0.788 |

![Overview of curated datasets]({{artifact:art_f6165da8-e2cb-44de-879e-aa8e7ecf2813}})

*(a) Size of the sets. (b) pIC50 distribution with the activity threshold
of 1 µM. (c) t-SNE of ECFP4 for the binding sets — each receptor occupies a
distinct chemical region. (d) Compounds tested across receptor pairs: the
overlap concentrates in ERα/ERβ (1,170 compounds) and THRα/THRβ (200), as
expected from the homology of the pairs.*

---

Each row of the CSV carries complete provenance: `molecule_chembl_id`, name,
original and standardized SMILES, InChIKey, kind/value/unit of measurement, pIC50,
pChEMBL, `assay_chembl_id` and textual description of the assay, `confidence_score`,
BAO format, organism, `target_chembl_id`, origin document, DOI, PMID, year and
journal. 83% of compounds have traceable DOI and PMID; the columns `all_assays`,
`all_docs`, `all_dois`, and `all_pmids` list **all** the sources that support the
aggregated value for each compound.

---
## 4. Modeling

**Representation.** Counted ECFP4 (radius 2, 2048 bits) concatenated with 201
RDKit 2D descriptors. For linear and margin-based models, the pipeline includes
`VarianceThreshold` and standardization — fitted only on the training set, inside
the cross-validation.

**Protocol.**

1. External split by **Bemis-Murcko scaffold** (80/20): the test set
   contains only molecular scaffolds absent from the training set. This is a
   deliberately severe test, close to the real-world scenario of prioritizing novel chemistry.
2. Model selection by grid search with **5-fold CV** within the training set
   (Ridge, SVR-RBF, Random Forest, XGBoost for regression; Logistic Regression,
   SVC-RBF, Random Forest, XGBoost for classification).
3. Refit on the full training set and evaluation on the external test set.
4. **y-scrambling** (5–10 permutations) on the best model.
5. **Applicability domain**: highest Tanimoto similarity (ECFP4) to the
   training set; the threshold is the 5th percentile of the nearest-neighbor
   similarity distribution within the training set.
6. **Production models** refitted on 100% of the curated data. The reported
   metrics continue to come from the protocol with the external test set.

**Classification.** Activity threshold at pIC50 ≥ 6 (IC50 ≤ 1 µM), a
common convention in endocrine-disruptor screening. The fraction of actives ranges from 0.41
(THRα) to 0.79 (functional AR).

| receptor   | assay     |   n_train_final | regression_model   |   Q2_cv |   R2_test_scaffold | classification_model   |   AD_threshold_tanimoto |
|:-----------|:----------|-----------------:|:-------------------|--------:|--------------------:|:-----------------------|---------------------:|
| AR         | binding   |             1807 | XGBoost            |   0.604 |               0.273 | XGBoost                |                0.5   |
| AR         | functional |              618 | RandomForest       |   0.45  |               0.076 | RandomForest           |                0.537 |
| ERα        | binding   |             2700 | XGBoost            |   0.793 |               0.66  | XGBoost                |                0.327 |
| ERα        | functional |              156 | RandomForest       |   0.589 |              -0.257 | RandomForest           |                0.324 |
| ERβ        | binding   |             1548 | XGBoost            |   0.602 |               0.098 | RandomForest           |                0.297 |
| PPARγ      | binding   |             1437 | XGBoost            |   0.662 |               0.477 | Logistic               |                0.555 |
| PPARγ      | functional |              228 | Ridge              |   0.627 |               0.376 | SVC                    |                0.494 |
| THRα       | binding   |              220 | RandomForest       |   0.72  |              -0.065 | RandomForest           |                0.571 |
| THRβ       | binding   |              396 | XGBoost            |   0.767 |               0.58  | XGBoost                |                0.619 |

---
## 5. Performance

![Model validation]({{artifact:art_69bd504f-ddad-4ae1-985f-cce7a0ed11e0}})

### pIC50 Regression

| receptor   | assay     | model_name   |   Q2_cv |   RMSE_cv |   R2_test |   RMSE_test |   MAE_test |
|:-----------|:----------|:-------------|--------:|----------:|-----------:|-------------:|------------:|
| AR         | binding   | XGBoost      |   0.604 |     0.632 |      0.273 |        0.776 |       0.604 |
| AR         | functional | RandomForest |   0.45  |     0.634 |      0.076 |        0.662 |       0.523 |
| ERα        | binding   | XGBoost      |   0.793 |     0.717 |      0.66  |        0.747 |       0.554 |
| ERα        | functional | RandomForest |   0.589 |     0.759 |     -0.257 |        1.068 |       0.867 |
| ERβ        | binding   | XGBoost      |   0.602 |     0.771 |      0.098 |        0.961 |       0.779 |
| PPARγ      | binding   | XGBoost      |   0.662 |     0.671 |      0.477 |        0.777 |       0.623 |
| PPARγ      | functional | Ridge        |   0.627 |     0.625 |      0.376 |        0.947 |       0.756 |
| THRα       | binding   | RandomForest |   0.72  |     0.535 |     -0.065 |        1.298 |       0.91  |
| THRβ       | binding   | XGBoost      |   0.767 |     0.697 |      0.58  |        0.427 |       0.339 |
### Active/inactive classification (pIC50 ≥ 6)

| receptor   | assay     | model_name   |   ROC_AUC_cv |   ROC_AUC_test |   PR_AUC_test |   BalAcc_test |   MCC_test |
|:-----------|:----------|:-------------|-------------:|----------------:|---------------:|---------------:|------------:|
| AR         | binding   | XGBoost      |        0.88  |           0.759 |          0.851 |          0.693 |       0.391 |
| AR         | functional | RandomForest |        0.844 |           0.772 |          0.984 |          0.5   |       0     |
| ERα        | binding   | XGBoost      |        0.951 |           0.855 |          0.965 |          0.749 |       0.485 |
| ERα        | functional | RandomForest |        0.872 |           0.94  |          0.974 |          0.875 |       0.831 |
| ERβ        | binding   | RandomForest |        0.894 |           0.821 |          0.965 |          0.748 |       0.391 |
| PPARγ      | binding   | Logistic     |        0.904 |           0.899 |          0.949 |          0.815 |       0.706 |
| PPARγ      | functional | SVC          |        0.909 |           0.932 |          0.949 |          0.91  |       0.821 |
| THRα       | binding   | RandomForest |        0.957 |           0.528 |          0.811 |          0.585 |       0.284 |
| THRβ       | binding   | XGBoost      |        0.918 |           0.841 |          0.736 |          0.783 |       0.486 |

### Independent external validation — BindingDB

817 compounds present in BindingDB and **absent** from ChEMBL, i.e., never
seen by the models in any step:

| receptor   |   n_external |   frac_in_domain |     R2 |   RMSE |   ROC_AUC |   BalAcc |   MCC |
|:-----------|------------:|------------------:|-------:|-------:|----------:|---------:|------:|
| ERα        |         378 |             0.944 |  0.502 |  1.272 |     0.883 |    0.703 | 0.56  |
| ERβ        |         235 |             0.702 | -0.286 |  1.286 |     0.761 |    0.77  | 0.519 |
| AR         |         128 |             0.68  |  0.466 |  0.917 |     0.837 |    0.811 | 0.682 |
| PPARγ      |          44 |             0.682 |  0.148 |  0.799 |     0.841 |    0.736 | 0.388 |
| THRα       |          17 |             0.824 |  0.888 |  0.359 |     1     |    1     | 1     |
| THRβ       |          15 |             0.133 |  0.362 |  0.689 |     1     |    0.964 | 0.681 |
### Reading the results

- **Y-scrambling rules out correlation by chance** across all endpoints: Q²
  between −0.43 and −0.14 with permuted labels, versus 0.45 to 0.79 with the
  real labels.
- **Classification generalizes better than regression.** Outside the training
  set, ROC-AUC is 0.76–0.94 in 8 of 9 endpoints, while regression R² drops
  from 0.45–0.79 (CV) to −0.26 to 0.66 (scaffold). For prioritizing novel
  compounds, use the **activity probability**; treat the predicted pIC50 as
  an approximate ranking, not as a quantitative estimate.
- **Four models are solid**: ERα binding (R² 0.66 on the scaffold test,
  RMSE 0.75 log), THRβ binding (0.58 / 0.43), PPARγ binding (0.48 / 0.78) and AR
  binding (0.27 / 0.78 — weak in regression, ROC-AUC 0.76 in classification).
- **Three models are fragile and should be used only as a qualitative alert**:
  ERα functional (n = 156), THRα binding (n = 220, ROC-AUC 0.53 on the
  scaffold test) and AR functional (balanced accuracy 0.50 — the classifier
  predicts the majority class, which represents 79% of the set).
- **THRα and THRβ have too much data in common and too little overall.** 200 of
  220 THRα compounds are also in THRβ; the models for the two receptors are not
  independent enough to support α/β selectivity prediction.
- **The applicability domain separates error clearly only for THRα**
  (RMSE 1.01 inside versus 1.98 outside). For the other endpoints the 5th
  percentile threshold is too permissive to work as a confidence filter —
  interpret the reported Tanimoto similarity as a continuous value, not as
  binary.

---
## 6. Using the predictor

```bash
python predict_edc.py --smiles "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1" --model-dir models_final
python predict_edc.py --input molecules.csv --smiles-col smiles --name-col name \
                      --model-dir models_final --out predictions.csv
```

For each molecule and each endpoint the output provides `pIC50_prev`, `prob_active`,
`tanimoto_nn` (similarity to the nearest neighbor in the training set), `in_domain`, and
a multi-target summary (`n_endpoints_active_in_domain`, `multitarget_profile`).

**Demonstration with known endocrine disruptors** (file
`predictions_demo.csv`) illustrates both the scope and the limitation of the model. The
bisphenol A is recovered as active in ERα and ERβ (predicted pIC50 5.96 and 6.02,
Tanimoto 1.0 — it is in the training sets) and flutamide as active in AR.
Genistein, however, is underestimated in ERα (4.94 predicted vs. ≈7.5 in the literature) and
rosiglitazone in PPARγ (4.75 vs. ≈8), **because neither of the two is in the
sets**: their potency is reported as EC50, Ki, or relative affinity, and
the IC50 cutoff excludes them. DDT, atrazine, and DEHP fall outside the applicability
domain in nearly all endpoints.

---
## 7. Limitations

1. **The IC50 cutoff is pharmacological, not toxicological.** ChEMBL's chemistry is
   dominated by drug-discovery series. Classic environmental contaminants
   (phthalates, organochlorines, triazines) are poorly represented and
   fall outside the applicability domain. For environmental screening, the natural
   path is to complement with the qHTS assays from Tox21/ToxCast — which report
   **AC50**, a distinct endpoint, and would therefore require a separate model. This expansion is detailed in section 9.
2. **Compounds missing by endpoint.** Well-known ligands are left out
   when their potency was published as EC50/Ki instead of IC50 (genistein,
   rosiglitazone). A parallel Ki dataset and another EC50 dataset would expand the
   coverage — without ever merging them with the IC50 one.
3. **Inherited experimental variability.** Even after removing discordant
   replicates, inter-laboratory dispersion imposes an error floor on the order
   of a few tenths of a log unit.
4. **Agonism and antagonism are not distinguished in the prediction.** The mode of action
   is annotated in the dataset, but the models predict potency, not the direction of the
   effect — a relevant point, since the toxicological consequence of agonizing or
   antagonizing a nuclear receptor is opposite.
5. **Absence of Cortellis** and of literature sources not indexed in ChEMBL
   or BindingDB.
6. **Intended use.** Prioritization of compounds for experimental assay. It does not
   replace in vitro testing nor constitute a regulatory risk assessment.

---

## 8. Reproducibility

| File | Content |
|---|---|
| `extract_chembl.py` | SQL extraction from the ChEMBL 37 SQLite dump |
| `extract_bindingdb.py` | filter of the complete BindingDB dump by accession |
| `curation.py` | RDKit standardization, conversion to pIC50, replicate aggregation |
| `pipeline_curate.py` | curation pipeline and dataset assembly |
| `features.py` | ECFP4 + descriptors, scaffold split, Tanimoto |
| `modeling.py` | training, selection, validation, y-scrambling, applicability domain |
| `finalize_models.py` | refit of production models on 100% of the data |
| `predict_edc.py` | multi-target prediction for new SMILES |

Environment: Python 3.11, RDKit 2026.03.6, scikit-learn 1.9.1, XGBoost, pandas,
numpy. Fixed seed (42) across all splits and models.

---

## 9. Scope expansion: independent model for AC50 (Tox21/ToxCast)

This section addresses limitation 7.1 and defines the scope of a new modeling effort. It is a
proposal: no results below have been obtained yet.
### 9.1 Motivation

The ChEMBL IC50 subset is dominated by drug-discovery series, and environmental contaminants
(phthalates, organochlorines, triazines) fall outside the applicability domain. The large-scale
quantitative high-throughput screening (qHTS) assays from Tox21 and ToxCast cover precisely
environmental and industrial chemistry. They measure **AC50**, the concentration at which the
response reaches 50% of its maximum amplitude.

### 9.2 Current status in the project

The existing Tox21 layer (`pipeline_tox21.py`, `modeling_tox21.py`, `models_tox21/`) contains
active/inactive classifiers per assay. AC50 is converted to pAC50 (6 − log10 AC50 in µM) and is only
used when the outcome is "Active". There is not yet a model that **predicts the value** of pAC50. ToxCast
(invitrodb) has also not been incorporated. The proposed expansion covers these two gaps.

### 9.3 Why AC50 cannot be merged into the IC50 model

1. **Nature of the parameter.** ChEMBL IC50 describes inhibition in curated literature assays, binding
   or functional, with heterogeneous protocols. AC50 is the parameter of a curve fit (Hill) over a fixed
   concentration grid, generated by a uniform screening pipeline. The two numbers are not on the same
   scale of meaning, even when the units coincide.
2. **Direction of the effect.** qHTS AC50 is reported for activation (agonist mode) and for inhibition
   (antagonist mode). The same AC50 can mean opposite effects, and the IC50 model does not distinguish
   direction (limitation 7.4).
3. **Data censoring.** Inactive substances do not have a well-defined AC50, or have an artificial limit
   value. The observed pAC50 distribution is truncated by the tested concentration range.
4. **Confounding by cytotoxicity.** In antagonist mode, a drop in signal may reflect cell death. Tox21
   provides viability counter-assays, which the IC50 layer does not have.
5. **Technology and cell line.** Luciferase or beta-lactamase reporter, cell line, and presence of
   autofluorescence introduce assay-specific bias, absent from binding data.
6. **Pooling effect.** Mixing the two regression targets would shift the target distribution,
   contaminate the applicability domain, and invalidate the conformal calibration, which is done per
   layer. The project already follows the rule of never merging IC50, Ki, Kd, and EC50, and AC50 follows
   the same rule.
### 9.4 Proposed architecture (independent of the IC50 layer)

| Component | Proposal |
|---|---|
| Modeling unit | one model per assay and mode (agonist or antagonist), without pooling across assays |
| Two-stage structure | stage 1: active/inactive classification (already exists, will be reused as a starting point). Stage 2: pAC50 regression only on the actives |
| Censoring | treat limit values as censored data (for example, regression with censoring or explicit exclusion with reporting), never as an observed value |
| Quality filters | exclude actives confounded by cytotoxicity or autofluorescence, using the counter-assays already present in `pipeline_tox21.py` |
| Features | the same as `features.py` (count-based ECFP4 and RDKit descriptors), to maintain operational compatibility |
| Validation | external Bemis-Murcko scaffold split (80/20), 5-fold CV in selection, y-scrambling, as in the other layers |
| Uncertainty and domain | conformal prediction and applicability domain calibrated **only** on the AC50 data |
| Output | new layer `models_ac50/`, with its own column in `_detail.csv`; never enters the average of the active fraction of the IC50 layers without an explicit rule |

### 9.5 Data

- **Tox21:** already used via PubChem, with the assays listed in `datasets_tox21_by_assay/`.
- **ToxCast (invitrodb):** to be incorporated. The exact version and the mapping of assay components to
  the receptors in scope (AR, ESR1, ESR2, PPARG, THRA, THRB) must be confirmed at the source before
  download. Comparability between assays will be verified using the same criterion as the other layers:
  same target, same effect direction, same technology, and nothing merged for convenience.
- **Curation:** the same standardization as `curation.py` (salts and solvents, neutralization, tautomer,
  InChIKey), with deduplication by InChIKey and provenance recording per measurement.

### 9.6 Acceptance criteria

1. Performance on the external scaffold-based test reported per assay, with interval and comparison against
   a trivial reference model (training mean).
2. Domain coverage measured on an environmental panel (phthalates, organochlorines, triazines, bisphenols),
   compared with the current coverage of the IC50 layers. Hypothesis to be tested: the coverage gain
   occurs in these classes.
3. Positive controls: the panel's positive controls (rosiglitazone in PPARG, genistein in
   ESR1) must continue to show a high alert after the integration.
4. No-merge test: verify that the new layer does not change the existing alert summary when
   switched off.
### 9.7 Risks

- Few actives per assay may render the regression of some receptor-mode pairs unfeasible.
- Tox21 and ToxCast assays for the same target may diverge, which requires analyzing the agreement
  before any union between the two sources.
- The AC50 of cell activation is not equivalent to binding affinity, and this must be made explicit in the
  interpretation of the alert.
