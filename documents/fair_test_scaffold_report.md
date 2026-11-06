# Fair test: the effect of putting large scaffold groups in training

Working copy in the sandbox; nothing in `~/Downloads/pipeline_edc_qsar` was changed. The current pipeline remains as it is.

## 1. Summary

Question: is it possible to improve the metrics by including the large scaffold groups in training? I answered this by measuring the three training runs on the **same test set**, which the previous re-evaluation did not do (each split rule generated a different test set, and the metrics were not comparable).

1. **There is a small gain, and it comes from the amount of data, not from the large groups themselves.** With the large groups in training (condition B, 80% of the data), the median RMSE drops by 0.018 (95% CI: 0.003 to 0.030; 19 of 26 endpoints improve) and the R² rises by 0.021. Relative to the median RMSE of 0.87, this is about 2%.
2. **Composition alone changes nothing.** When the training set with large groups is reduced to the size of the training set without them (condition C), the differences disappear: RMSE −0.007 (95% CI: −0.019 to 0.016), R² −0.003 (95% CI: −0.018 to 0.012).
3. **Classification:** the ROC-AUC of the ChEMBL endpoints rises by 0.007 with B (95% CI: −0.001 to 0.015), and that of Tox21 by 0.006 (95% CI: −0.003 to 0.018). Both intervals include zero.
4. **The production model does not change.** It is already trained with 100% of the data, large groups included. The split rule only changes how performance is estimated.
5. **The R² jump from the previous re-evaluation (median 0.02 to 0.41) is not a model gain.** On the same test set, the R² gain is 0.02. The remainder is consistent with the composition of the old test set, which had a restricted activity range (measured in the previous re-evaluation); I did not redo that decomposition here.
## 2. Design

- **Large groups (L):** the scaffold groups that `features.scaffold_split` sends to the test set (greedy, largest to smallest, up to 20% of the compounds). **S** is the remainder.
- **Common test set (T):** about 20% of the compounds, sampled by scaffold group within S, with a random seed. No compound in T appears in any training set. Median per data group: IC50: 124 compounds from 50 scaffolds; Ki/Kd/EC50: 54 compounds from 21 scaffolds; Tox21: 991 compounds from 750 scaffolds. **R** = S minus T.
- **Training conditions**, all assessed on T:
  - **A:** training = R (large groups outside; 60% of the data).
  - **B:** training = R ∪ L (large groups inside; 80% of the data).
  - **C:** subsample of R ∪ L with the size of R (isolates composition from amount of data).
- **Models:** the algorithm and hyperparameters that the protocol selected under the old regime, fixed across the three conditions (I did not rerun the CV-based selection). 26 ChEMBL endpoints (9 IC50 and 17 Ki/Kd/EC50; 24 with a classifier) and 11 Tox21 endpoints (classification only).
- **Test draws:** 4 for ChEMBL (seeds 11 to 14) and 3 for Tox21 (seeds 11 to 13). For each endpoint, I use the median across draws.
- **Criteria fixed before looking at the results:** unit of analysis = endpoint; paired Wilcoxon test; 95% CI by bootstrap across endpoints; exclusion of any endpoint-draw pair with R² < −5 (numerical failure). No row was excluded by this criterion.

| Group | Condition | Training n (median) | Fraction of training coming from large groups | Maximum Tanimoto of test to training (median) |
|---|---|---|---|---|
| IC50 | A | 370 | 0.00 | 0.652 |
| IC50 | B | 494 | 0.25 | 0.659 |
| IC50 | C | 370 | 0.25 | 0.643 |
| Ki/Kd/EC50 | A | 163 | 0.00 | 0.635 |
| Ki/Kd/EC50 | B | 217 | 0.25 | 0.644 |
| Ki/Kd/EC50 | C | 163 | 0.25 | 0.635 |
| Tox21 | A | 2975 | 0.00 | 0.378 |
| Tox21 | B | 3966 | 0.25 | 0.404 |
| Tox21 | C | 2975 | 0.25 | 0.386 |
## 3. Results

Δ = condition − reference. For RMSE, a negative Δ is better; for R² and ROC-AUC, positive is better. The count columns show how many endpoints have a positive and negative difference.

| Metric | Comparison | Endpoints | Median Δ | 95% CI | Differences > 0 / < 0 | p (Wilcoxon) |
|---|---|---|---|---|---|---|
| RMSE (regression) | B - A | 26 | −0.018 | −0.030 to −0.003 | 7 / 19 | 0.007 |
| RMSE (regression) | C - A | 26 | −0.007 | −0.019 to 0.016 | 12 / 14 | 0.727 |
| RMSE (regression) | B - C | 26 | −0.019 | −0.038 to −0.001 | 8 / 18 | 0.056 |
| R² (regression) | B - A | 26 | 0.021 | 0.002 to 0.036 | 20 / 6 | 0.024 |
| R² (regression) | C - A | 26 | −0.003 | −0.018 to 0.012 | 12 / 14 | 0.784 |
| R² (regression) | B - C | 26 | 0.027 | 0.016 to 0.039 | 20 / 6 | 0.038 |
| ROC-AUC (ChEMBL) | B - A | 24 | 0.007 | −0.001 to 0.015 | 16 / 8 | 0.031 |
| ROC-AUC (ChEMBL) | C - A | 24 | −0.005 | −0.012 to 0.004 | 10 / 14 | 0.422 |
| ROC-AUC (ChEMBL) | B - C | 24 | 0.012 | 0.008 to 0.022 | 18 / 6 | 0.034 |
| ROC-AUC (Tox21) | B - A | 11 | 0.006 | −0.003 to 0.018 | 6 / 5 | 0.123 |
| ROC-AUC (Tox21) | C - A | 11 | −0.010 | −0.018 to 0.002 | 3 / 8 | 0.365 |
| ROC-AUC (Tox21) | B - C | 11 | 0.008 | −0.002 to 0.020 | 8 / 3 | 0.019 |

- **Consistency across draws:** condition B has lower RMSE than A in 66% of the 104 endpoint-draw pairs, and C in 45%. For ROC-AUC, B is better than A in 62% of the ChEMBL pairs and in 70% of the Tox21 pairs (C: 45% and 39%).
- **Endpoint size:** the gain of B over A does not correlate detectably with training set size (Spearman −0.16 for RMSE and −0.22 for ROC-AUC, p > 0.3; n = 26 and 24).
- **The test set gets only slightly closer to the training set:** the median maximum Tanimoto similarity of test to training, for example, goes from 0.652 to 0.659 for IC50 and from 0.378 to 0.404 for Tox21.

![Gain per endpoint, conditions B and C relative to A]({{artifact:art_ea568f8f-7232-4974-8e25-fbe9fe645b42}})
## 4. Interpretation

1. Adding the large groups to the training set improves performance on compounds with novel scaffolds, but only slightly: about 0.02 pAct unit in RMSE, 0.02 in R², and 0.007 in ROC-AUC across the ChEMBL endpoints.
2. The gain is explained by the increase from 60% to 80% of the data: B − C (quantity only) has a Δ of −0.019 in RMSE, 0.027 in R², and 0.012 in ROC-AUC, and C − A (composition only) is close to zero. There is no evidence that the large groups are more informative per compound than the small ones, although the CI of the C − A contrast still allows for an effect of up to about 0.02 in RMSE.
3. Regarding the decision about the `scaffold_split` rule: adopting large groups in training is defensible as an evaluation protocol (training on 80% and testing on less populous scaffolds, similar to prioritizing novel compounds), but it is not a lever for improving the model. Larger gains would require other approaches, which I did not test (for example, model selection via scaffold-grouped CV and transformation of the `Ipc` descriptor).
4. When comparing metrics between split rules, compare on a common test set. The R² figures for the corrected rule and the current rule are not directly comparable.
## 5. Limitations

1. **Test only with small scaffolds.** T is sampled within S, so the conclusion holds for new, sparsely populated scaffolds, and not for the dominant chemical families (L). The test of the current rule (extrapolation to L) is not comparable.
2. **Algorithm and hyperparameter selection used CV on S, which contains T.** This may slightly inflate the absolute metrics of the three conditions. The differences between conditions are less affected, because the configuration is the same, but the effect was not quantified. The configuration was not re-optimized for condition B.
3. **Sizes different from those used before.** The training set of condition A has 60% of the data, so the absolute metrics are not comparable with those of the previous report (80%).
4. **Few random draws** (4 on ChEMBL, 3 on Tox21), and the fifth seed planned for ChEMBL did not run due to an argument list error in the queue.
5. **Non-independent endpoints.** IC50, Ki, Kd and EC50 for the same receptor share compounds, so the CIs between endpoints tend to be optimistic.
6. **No correction for multiple comparisons** (24 tests in total: 8 metrics and 3 contrasts). The isolated p-values are indicative only. Only the B − A and B − C contrasts point in the same direction across several metrics.
7. **Tox21:** the CSVs in the folder do not match the deployed models (divergence already recorded); this test applies to the CSVs.
8. **Ridge failures with `Ipc`:** no numerical failure appeared here, although 3 regressions use Ridge in the old regime. I did not investigate why the problem, previously seen with the corrected rule, does not recur in this test; this does not eliminate the risk.

## 6. Files

`fair_test.py` runs the test per endpoint and `fair_test_analysis.py` generates the tables; `fair_long.csv` provides one row per endpoint, condition, task and random draw; `fair_summary_paired.csv` provides all comparisons, including MAE, Spearman and PR-AUC.
### All comparisons

| Metric | Comparison | Endpoints | Median Δ | 95% CI | Differences > 0 / < 0 | p (Wilcoxon) |
|---|---|---|---|---|---|---|
| RMSE (regression) | B - A | 26 | −0.018 | −0.030 to −0.003 | 7 / 19 | 0.007 |
| RMSE (regression) | C - A | 26 | −0.007 | −0.019 to 0.016 | 12 / 14 | 0.727 |
| RMSE (regression) | B - C | 26 | −0.019 | −0.038 to −0.001 | 8 / 18 | 0.056 |
| MAE (regression) | B - A | 26 | −0.012 | −0.023 to 0.007 | 10 / 16 | 0.075 |
| MAE (regression) | C - A | 26 | −0.000 | −0.012 to 0.022 | 13 / 13 | 0.764 |
| MAE (regression) | B - C | 26 | −0.016 | −0.031 to 0.016 | 10 / 16 | 0.129 |
| R² (regression) | B - A | 26 | 0.021 | 0.002 to 0.036 | 20 / 6 | 0.024 |
| R² (regression) | C - A | 26 | −0.003 | −0.018 to 0.012 | 12 / 14 | 0.784 |
| R² (regression) | B - C | 26 | 0.027 | 0.016 to 0.039 | 20 / 6 | 0.038 |
| Spearman (regression) | B - A | 26 | 0.018 | 0.000 to 0.026 | 18 / 8 | 0.080 |
| Spearman (regression) | C - A | 26 | −0.008 | −0.018 to 0.003 | 9 / 17 | 0.280 |
| Spearman (regression) | B - C | 26 | 0.019 | 0.010 to 0.029 | 20 / 6 | 0.025 |
| ROC-AUC (ChEMBL) | B - A | 24 | 0.007 | −0.001 to 0.015 | 16 / 8 | 0.031 |
| ROC-AUC (ChEMBL) | C - A | 24 | −0.005 | −0.012 to 0.004 | 10 / 14 | 0.422 |
| ROC-AUC (ChEMBL) | B - C | 24 | 0.012 | 0.008 to 0.022 | 18 / 6 | 0.034 |
| ROC-AUC (Tox21) | B - A | 11 | 0.006 | −0.003 to 0.018 | 6 / 5 | 0.123 |
| ROC-AUC (Tox21) | C - A | 11 | −0.010 | −0.018 to 0.002 | 3 / 8 | 0.365 |
| ROC-AUC (Tox21) | B - C | 11 | 0.008 | −0.002 to 0.020 | 8 / 3 | 0.019 |
| PR-AUC (ChEMBL) | B - A | 24 | 0.007 | 0.003 to 0.014 | 20 / 4 | 0.001 |
| PR-AUC (ChEMBL) | C - A | 24 | −0.000 | −0.006 to 0.009 | 12 / 12 | 0.546 |
| PR-AUC (ChEMBL) | B - C | 24 | 0.007 | 0.003 to 0.011 | 19 / 5 | 0.016 |
| PR-AUC (Tox21) | B - A | 11 | 0.012 | −0.004 to 0.030 | 8 / 3 | 0.032 |
| PR-AUC (Tox21) | C - A | 11 | 0.002 | −0.025 to 0.023 | 6 / 5 | 0.966 |
| PR-AUC (Tox21) | B - C | 11 | 0.011 | 0.006 to 0.024 | 10 / 1 | 0.002 |
