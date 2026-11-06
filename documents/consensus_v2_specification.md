# Consensus v2: two-axis report

Implementation of step 1 of the proposed architecture: **aggregation-only changes**. No model was
retrained or altered, no individual vote changed, and nothing in `~/Downloads/pipeline_edc_qsar` was
touched. The `consensus_v2.py` reads the `<output>_detail.csv` that your `predict_edc2.py` already produces and
rewrites only how these votes are summarized.

## 1. What changed

| | `_summary.csv` current (v1) | Consensus v2 |
|---|---|---|
| Evidence unit | the model (regression and classification for the same endpoint count twice) | the endpoint (1 piece of evidence) |
| Role of each task | both vote equally | classification decides positivity; regression gives potency with a conformal interval |
| Aggregation | a single level, summing ChEMBL and Tox21 | two axes, never summed: pharmacological and screening |
| Within the axis | fraction of active votes | best evidence and how many endpoints corroborate it |
| Model outside the domain | disappears and dilutes the fraction | explicit `not_assessed` state, with coverage alongside |
| Direction of effect | absent | label derived from the mode of the Tox21 assays |
| Pairs without coverage | no row in the file | row with `no coverage` |

Why classification decides: outside training it generalizes better (test ROC-AUC from 0.76 to 0.94
on the IC50 endpoints, versus test R² from −0.26 to 0.66). Regression remains in the report as
an order of magnitude for potency, with the conformal interval, and not as a second vote.

Why best evidence instead of the fraction: the endpoints are not independent. Ki and EC50 share
compounds (205 of the 271 in `AR_EC50_F` are also in `AR_Ki_B`), and the 4 Tox21 assays for AR and ERα are
2 modes in 2 cell lines. The fraction gives implicit weight to the receptor that has more datasets.

## 2. States and reading

Each axis receives a state: `positive`, `negative`, `not_assessed` (a deployable model exists, but the
molecule is outside the domain), or `no_model` (no deployable model exists for that receptor
on that axis). The combined reading names the 9 possible combinations, without summing anything.
## 3. Coverage by receptor

This is the cap for each receptor, the same for any molecule. It explains why a
negative result in ERβ or in THRβ is worth less than in AR or ERα.

| Receptor | Pharmacological endpoints | Tox21 assays | Measures | Tox21 modes |
|---|---|---|---|---|
| AR | 3 | 4 | EC50, IC50, Ki | agonist, antagonist |
| ESR1 | 4 | 4 | EC50, IC50, Ki | agonist, antagonist |
| ESR2 | 4 | 0 | EC50, IC50, Ki | — |
| PPARG | 5 | 2 | EC50, IC50, Ki | agonist, antagonist |
| THRB | 1 | 0 | IC50 | — |
| TR | 0 | 1 | — | antagonist |

THRα does not appear: its 6 ChEMBL models are all non-deployable and there is no Tox21 assay with
that label. ERβ has no Tox21 assay, and the Tox21 `TR` label does not resolve α from β.
## 4. Comparison with the current summary, on the reference panel

I ran your `predict_edc2.py` without changes on `panel_ref.csv` (18 molecules) and re-aggregated the
`_detail.csv`. v1 produced 83 rows, and v2 produces 126 molecule × receptor pairs
(43 of them had no row in v1).

| v2 reading | high | moderate | low | no alert | (no row in v1) |
|---|---|---|---|---|---|
| convergent negative | 0 | 0 | 0 | 1 | 0 |
| convergent positive | 10 | 0 | 0 | 0 | 0 |
| divergent: pharmacological positive, screening negative | 0 | 0 | 1 | 0 | 0 |
| divergent: screening positive, pharmacological negative | 0 | 1 | 3 | 0 | 0 |
| negative only on the screening axis (partial coverage) | 0 | 0 | 0 | 34 | 0 |
| negative only on the pharmacological axis (partial coverage) | 0 | 0 | 0 | 5 | 0 |
| no coverage | 0 | 0 | 0 | 0 | 43 |
| pharmacological positive only (no screening coverage) | 0 | 4 | 2 | 0 | 0 |
| screening positive only (no pharmacological coverage) | 0 | 5 | 17 | 0 | 0 |

What the table shows:

1. **The 10 "high" pairs in v1 all turn into "convergent positive".** The strongest level in v1 already was, in fact, convergence between the two axes.
2. **The 10 "moderate" pairs split into three distinct situations**, which the single level lumped together:
   5 are screening positive only with no pharmacological coverage at all (Methoxychlor, Vinclozolin, and
   Triclosan in AR; Triclosan and Flutamide in ERα), 4 are pharmacological positive only with no screening
   coverage (the four for ERβ, which has no Tox21 assay), and 1 is a divergence (Zearalenone in ERα, with
   screening positive and pharmacological negative).
3. **Of the 40 "no alert" pairs in v1, only 1 is negative on both axes.** The other 39 are negative on
   one axis with the other lacking coverage. In practice, "no alert" almost never meant that both
   lines of evidence said no.
4. **The coverage bias becomes visible instead of embedded.** The 5 moderates for AR and ERα with no
   pharmacological coverage were indistinguishable, in v1, from a moderate supported by both layers.
## 5. Limitations

1. **This does not improve accuracy, and it was not validated for accuracy.** It is a reporting change: it makes auditable what the single level was hiding. The panel has 15 positives and 3 negatives assessed, with the class assigned per molecule and not per receptor, so it is not suitable for comparing v1 and v2 by sensitivity. For that, external labels per receptor are missing, which the screening axis does not have (the Tox21 label is itself the training signal).
2. **The probability thresholds remain those of each model** (0.11 to 0.40 in Tox21 versus 0.5 in ChEMBL). Calibration is step 2 of the proposal and is not implemented; until then, the probability of one axis is not comparable with that of the other.
3. **Direction comes only from Tox21.** The pharmacological axis remains as `undetermined`, because the mode annotation in ChEMBL is inferred by regular expression and is not reliable in binding assays.
4. **v2 does not produce a single level and does not rank molecules.** If you need a prioritization queue, that ordering rule is a separate decision, which I did not make here.
5. **The panel consists of environmental chemicals**, almost all outside the ChEMBL domain. The proportion of `not_assessed` in the pharmacological axis would be lower with drug-like molecules. I did not measure this.
6. **Potency is not a single scale.** pIC50, pKi, and pEC50 are not interchangeable. In 15 of the 27 pairs with potency in the pharmacological axis there is more than one measure available, so `pAct_best` is the maximum across different scales and only makes sense read together with `pAct_best_measure`. Use the `potency_by_measure` column, which provides the potency of each measure separately.
7. **In the 504 deployable endpoints of this panel, the decision always came from classification
   when there was coverage** (257 cases) or from no
   task (247 outside the domain). The alternative path
   (endpoint with deployable regression and without classifier) exists in the code, but was not exercised
   here, so it is not tested.

## 6. How to use

```
# 1. prediction with the original pipeline, from inside the pipeline folder (the default paths matter:
#    predict_edc2.py identifies the layer by the folder name, so absolute --model-dirs breaks the label)
cd ~/Downloads/pipeline_edc_qsar
python predict_edc2.py --input panel_ref.csv --name-col name --out /path/to/output/panel
# 2. reaggregation along two axes
python consensus_v2.py /caminho/output/panel_detail.csv /caminho/output/consensus_v2
```

Outputs: `_endpoints.csv` (one row per molecule × receptor × endpoint), `_summary.csv` (one row per
molecule × receptor, with the two axes and the coverage) and `_coverage.csv` (the cap of each receptor).
