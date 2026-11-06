"""Paired analysis of the fair test (reads out/fair_*.json and writes out/fair_long.csv and out/fair_summary_paired.csv).

Criteria defined before looking at the results:
  - unit of analysis = endpoint; per endpoint, median across test samplings;
  - samplings used: 11-14 (ChEMBL) and 11-13 (Tox21);
  - regression: excludes the pair (endpoint, sampling) if R² < -5 in any condition (numerical failure);
  - paired comparison via Wilcoxon; 95% CI of the median of differences via bootstrap across endpoints (4,000 replicates, seed 0).
Usage: cd work; python fair_test_analysis.py
"""
import glob, json
import numpy as np, pandas as pd
from scipy import stats

C = ["A_no_large_groups", "B_with_large_groups", "C_equal_size"]
R = []
for f in sorted(glob.glob("out/fair_*.json")):
    R += json.load(open(f))
D = pd.DataFrame(R)
D = D[~((D.group != "tox21") & (D.random_seed == 15))]
D.to_csv("out/fair_long.csv", index=False)

reg = D[D.task == "regression"]
failure = reg.groupby(["group", "endpoint", "random_seed"]).r2.transform("min") < -5
print("regression rows with numerical failure:", int(failure.sum()), "of", len(reg))
regok, cl = reg[~failure], D[D.task == "classification"]
rng = np.random.default_rng(0)


# Paired comparison of one metric between training conditions: Wilcoxon plus a bootstrap CI of the median difference.
def paired(d, metric, groups):
    m = d[d.group.isin(groups)].groupby(["endpoint", "condition"])[metric].median().unstack().dropna(subset=C)
    out = []
    for a, b, name in (("A_no_large_groups", "B_with_large_groups", "B - A"), ("A_no_large_groups", "C_equal_size", "C - A"),
                       ("C_equal_size", "B_with_large_groups", "B - C")):
        dlt = m[b] - m[a]
        n = len(dlt)
        p = stats.wilcoxon(dlt).pvalue if n >= 6 and (dlt != 0).any() else np.nan
        bs = [np.median(rng.choice(dlt.values, n, replace=True)) for _ in range(4000)]
        out.append(dict(comparison=name, n_endpoints=n, median_ref=m[a].median(), median_alt=m[b].median(),
                        median_delta=dlt.median(), ic95_lower=np.percentile(bs, 2.5), ic95_upper=np.percentile(bs, 97.5),
                        n_positive_diff=int((dlt > 0).sum()), n_negative_diff=int((dlt < 0).sum()), p_wilcoxon=p))
    return pd.DataFrame(out)


tabs = []
for name, d, met, gr in [("RMSE (regression)", regok, "rmse", ["ic50", "layers"]), ("MAE (regression)", regok, "mae", ["ic50", "layers"]),
                         ("R² (regression)", regok, "r2", ["ic50", "layers"]), ("Spearman (regression)", regok, "spearman", ["ic50", "layers"]),
                         ("ROC-AUC (ChEMBL)", cl, "roc_auc", ["ic50", "layers"]), ("ROC-AUC (Tox21)", cl, "roc_auc", ["tox21"]),
                         ("PR-AUC (ChEMBL)", cl, "pr_auc", ["ic50", "layers"]), ("PR-AUC (Tox21)", cl, "pr_auc", ["tox21"])]:
    t = paired(d, met, gr)
    t.insert(0, "metric", name)
    tabs.append(t)
pd.concat(tabs).to_csv("out/fair_summary_paired.csv", index=False)


# Fraction of endpoint x seed pairs in which each condition beats condition A.
def consist(d, metric, lower_is_better):
    w = d.pivot_table(index=["endpoint", "random_seed"], columns="condition", values=metric).dropna()
    f = lambda a, b: (w[b] < w[a]) if lower_is_better else (w[b] > w[a])
    return dict(metric=metric, pairs=len(w), B_better_than_A=round(f(C[0], C[1]).mean(), 3), C_better_than_A=round(f(C[0], C[2]).mean(), 3))


pd.DataFrame([dict(group="regression", **consist(regok, "rmse", True)), dict(group="ChEMBL", **consist(cl[cl.group != "tox21"], "roc_auc", False)),
              dict(group="Tox21", **consist(cl[cl.group == "tox21"], "roc_auc", False))]).to_csv("out/fair_consistency_random_seeds.csv", index=False)
print("OK")
