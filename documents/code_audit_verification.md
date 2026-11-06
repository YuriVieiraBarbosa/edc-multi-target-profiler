# Verification of the code audit

I read the excerpts cited in `pipeline_edc_qsar` and ran diagnostics on the 38 endpoints with a dataset in the folder
(9 IC50, 17 Ki/Kd/EC50, and 12 Tox21). Nothing in your folder was altered. What follows separates what
I confirmed, what needs a caveat, and what the audit underestimated.

## 1. Items I confirmed in the code

| Item | Evidence |
|---|---|
| A1 | `extract_bindingdb.py:33` writes `cache/bindingdb_ic50_raw.parquet`; `pipeline_curate.py:85` reads `.pkl`. |
| A2 | `bindingdb_std_lookup.pkl` appears only in `pipeline_curate.py:86`, in a read. No script writes it. |
| A3 | `README_calibration.csv` has 22 columns; `conformal.py` writes `results/conformal_and_thresholds.csv`. The 5 columns (`production_threshold`, `threshold_gain`, `conformal_mode`, `coverage_scaffold`, `not_deployable`) come from a script that is not in the folder. |
| A5 | `data/` and `results/` do not exist. `finalize_models.py:12,15` and `conformal.py` depend on both. |
| A7 | `run_edc.py:28,48` call `subprocess` without `cwd`. |
| B1 | `predict_edc2.py:32` compares the folder name as text. Already documented in three places. |
| B7 | The `agg` uses `ev.loc[s.index, "active"]`, dependent on index alignment. |
| C1 | The docstring says "larger scaffolds go to training"; the loop fills the **test** set first. |
| C2–C5 | `conformal.py:50-57` calibrates on a model fitted on 75% of the training set; `finalize_models.py:26` fits on 100%. The domain threshold is recalculated over 100% (`finalize_models.py:27-28`), and `Q2_cv`/`R2_test_scaffold` from the 80% model are written inside the production artifact. |

**C2 closes a gap of mine.** In the fair-test report I measured that conformal undercoverage
persists under both split rules, and said I didn't know the cause. The likely cause is here: the quantile
is calibrated on one model and applied to another, fitted on more data. This hasn't been tested, but it is
an explanation that the split did not account for.

## 2. Item that needs a caveat

**A4 is not a code bug.** `pipeline_layers.py:106` writes to `data_layers/` and `modeling_layers.py:95`
reads from `data_layers/`; the same holds for `data_tox21/`. The code is internally consistent. The
discrepancy is that in the delivered folder these outputs appear as `datasets_layers_by_endpoint/` and
`datasets_tox21_by_assay/`, meaning they were renamed after being generated. The fix is a packaging
issue (one name, in one place), not a logic issue.
## 3. B4: the fallback exists, but it never triggers

`features.py:73-75` falls back to random split when the greedy test falls below half of the target, without
logging anything. I confirmed the code. I ran the greedy loop on the 38 endpoints:

- **the fallback does not trigger in any of them.** In all of them, greedy filling reaches exactly the
  target size of the test (ratio of 1.00 between achieved and target, in the worst case).

So no endpoint was randomly validated, and the fair `scaffold_split` test I ran was not
affected. This downgrades B4 from "endpoints may have been randomly validated" to a latent risk:
it's worth logging the flag, but there is no number to revise.

## 4. B5 is worse than the audit describes

Molecules with no Murcko scaffold (acyclic, or function failure) all fall into the `""` group. I measured the
size of this pseudo-group and where it lands in the split:

- **In Tox21 it is huge:** 18% to 23% of each assay, from 654 to 1,392 compounds. In 7 of the 12 assays it is
  the largest group in the dataset.
- **The effect is bimodal, and depends only on whether the group fits within the test budget or not.** In 15 of the 18
  endpoints with acyclics, the test ends up with **zero** acyclic molecules. In the other three, the
  opposite happens:

| Endpoint | Acyclic | test n | Fraction of test that is acyclic |
|---|---|---|---|
| `TOX21_AR_antagonist_HEK293T` | 944 | 992 | 0.95 |
| `TOX21_TR_antagonist_GH3` | 654 | 712 | 0.92 |
| `PPARG_Kd_B` | 4 | 13 | 0.31 |

In other words, in these two Tox21 assays the external test is, in practice, "all acyclic molecules", and
in the other nine it is "no acyclic molecule". The two situations are bad for opposite reasons, and the
difference between them is just the happenstance of `0.2 × n` being greater or smaller than the pseudo-group.

This directly affects the test ROC-AUC of these two endpoints, which were included in my comparisons. The
minimal fix is to treat each acyclic molecule as its own group, and that changes numbers.
## 5. Items in `consensus_v2.py` that are mine — fixed

- **E9** (order depended on the input): `groupby` now has `sort=True` and the pairs come out sorted.
- **E10** (missing columns): `load()` now creates `predicted_value`, `interval_lower`, `interval_upper`,
  `prob_active`, `threshold`, `tanimoto_nn` and `detail` when missing, so a `_detail.csv` without
  conformal calibration no longer breaks.
- **E11** (O(n²)): a `groupby` by `(name, receptor, axis)` replaced the re-filtering per pair.
- **E12** (silent `.iloc[0]`): now raises `ValueError` naming the endpoint and pointing to the usual
  cause, which is a repeated name in the input (the B2). This turns B2's silent-corruption path
  into an explicit error, at least in the reaggregation.

I verified that the content did not change: the 126 rows of the summary and the 666 of the endpoints file
have the same values as before, just in deterministic order. The `ValueError` was tested by duplicating a name.

**D3 is also partly mine.** The alert rule is copied in `panel_layer_ablation.py`. It exists
there to reproduce the `_summary.csv` and prove that the ablation starts from the same calculation, and this
reproduction is checked pair by pair. Even so, it is a third copy.

## 6. What I would prioritize

1. **A3.** As long as the post-processing script is not rebuilt, the 18 IC50 models cannot
   be regenerated by anything in the folder. It is the only item that prevents redoing what is deployed.
2. **B1 and B2.** These are the two paths to wrong results without warning. B2 now fails explicitly
   in the reaggregation, but it still merges rows in `predict_edc2.py`.
3. **B5.** It is not on the audit's priority list and, based on the numbers above, it should be: two
   Tox21 endpoints have a practically degenerate external test.
4. **C2 to C4.** These are decisions, but they explain two already-measured findings (conformal
   subcoverage and deployed domain wider than the validated one). It is worth freezing with the explanation recorded.

## 7. Caveats of this verification

- I read the cited excerpts, not the 24 entire files. Architecture and style items (D and E) were not
  audited one by one; I confirmed only the ones that touch `consensus_v2.py`.
- I did not test the proposed fixes for the **D** items, because they all change numbers and depend on your
  decision.
- The split diagnosis uses the datasets that are in the folder plus the project's IC50 artifact. If the
  datasets are regenerated with another RDKit version, tautomeric canonicalization changes the InChIKeys and
  the diagnosis needs to be redone (item C9, which I did not verify).
