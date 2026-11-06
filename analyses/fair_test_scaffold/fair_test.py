"""Fair test of the effect of putting large scaffold groups into training.

Problem with the previous comparison: each split rule produces a different TEST, so the metrics are not comparable.
Here the test T is the same across the three training conditions:

  L = large groups (the ones that features.scaffold_split sends to the test: greedily from largest to smallest, up to 20%).
  S = the rest (small groups and unique scaffolds).
  T = ~20% of the compounds, sampled PER scaffold GROUP within S (seed). R = S \ T.

  A  train = R                      (large groups outside training, as in the current rule's training composition)
  B  train = R ∪ L                  (large groups in training, as in the corrected rule)
  C  train = subsample of R ∪ L with the size of R (isolates composition from amount of data)

Algorithm and hyperparameters: those chosen by the protocol in the old regime (out/old_*.json), fixed across the three conditions.
Does not rerun CV selection. Usage: cd work; python fair_test.py <ic50|layers|tox21> <endpoint> <seed> [<seed> ...]
"""
import sys, os, json, time, warnings
from collections import defaultdict

warnings.filterwarnings("ignore")
sys.path.insert(0, os.getcwd())
import numpy as np, pandas as pd
from scipy import stats
from sklearn.base import clone
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error, roc_auc_score, average_precision_score

import features as F
import modeling as M
import modeling_tox21 as MT

NJ = int(os.environ.get("FAIR_NJOBS", "2"))
group, ep = sys.argv[1], sys.argv[2]
seeds = [int(s) for s in sys.argv[3:]]
FOLDER = {"ic50": ("data_ic50", "pIC50_final"), "layers": ("data_layers", "pAct_final"), "tox21": ("data_tox21", "specific_active")}[group]
out = f"out/fair_{group}_{ep}.json"
if os.path.exists(out):
    sys.exit(0)

df = pd.read_csv(f"{FOLDER[0]}/dataset_{ep}.csv")
cache = f"cache/{group}_{ep}.npz"
if os.path.exists(cache):
    z = np.load(cache, allow_pickle=True); X, bits, ok = z["X"], z["bits"], z["ok"]
else:
    X, bits, names, ok = F.featurize(df["std_smiles"].tolist())
    np.savez_compressed(cache, X=X, bits=bits, ok=ok)
df = df[ok].reset_index(drop=True)
sm = df["std_smiles"].tolist()
y = df[FOLDER[1]].to_numpy(dtype=float)
ycls = y.astype(int) if group == "tox21" else (y >= 6.0).astype(int)
N = len(df)
sc = np.array([F.murcko_scaffold(s) for s in sm])

S, L = F.scaffold_split(sm, 0.2, 42)          # S = training under the current rule; L = test under the current rule (large groups)
old = json.load(open(f"out/old_{group}_{ep}.json"))
res, rows = old["res"], old["rows"]


# Rebuild the algorithm and hyperparameters that the original protocol selected for this endpoint.
def config(task):
    if group == "ic50":
        name = res["best_regression"] if task == "regression" else res.get("best_classification")
        key = "best_params"
    elif group == "layers":
        name = res["best_model"] if task == "regression" else res.get("best_classifier")
        key = "params"
    else:
        name, key = res["best_model"], "params"
    if not name:
        return None
    r = [x for x in rows if x["model_name"] == name and (group == "tox21" or x["task"] == task)][0]
    zoo = (M.reg_zoo(X.shape[1]) if task == "regression" else (MT.zoo() if group == "tox21" else M.clf_zoo()))
    est = clone(zoo[name][0]).set_params(**json.loads(r[key]))
    if "n_jobs" in est.get_params():
        est.set_params(n_jobs=NJ)
    return name, est


tasks = ["classification"] if group == "tox21" else ["regression", "classification"]
cfg = {t: config(t) for t in tasks}
records = []
t0 = time.time()
for seed in seeds:
    rng = np.random.default_rng(seed)
    grp = defaultdict(list)
    for i in S:
        grp[sc[i]].append(i)
    ch = list(grp.values())
    target = int(round(0.2 * N)); T = []
    for k in rng.permutation(len(ch)):
        if len(T) + len(ch[k]) <= target:
            T.extend(ch[k])
        if len(T) >= target:
            break
    T = np.array(sorted(T)); Tset = set(T.tolist())
    R = np.array([i for i in S if i not in Tset]); R_plus_L = np.concatenate([R, L])
    C = np.sort(np.random.default_rng(seed + 1000).choice(R_plus_L, size=len(R), replace=False))
    conditions = {"A_no_large_groups": R, "B_with_large_groups": R_plus_L, "C_equal_size": C}
    for cn, tr in conditions.items():
        simT = F.tanimoto_max_to_train(bits[T], bits[tr])
        base = dict(group=group, endpoint=ep, random_seed=seed, condition=cn, n_train=int(len(tr)), n_test=int(len(T)),
                    train_frac_of_L=float(np.isin(tr, L).mean()), test_scaffolds=int(len(set(sc[T]))),
                    tanimoto_max_test_median=float(np.median(simT)), test_prevalence=float(ycls[T].mean()),
                    train_prevalence=float(ycls[tr].mean()))
        for task in tasks:
            if cfg[task] is None:
                continue
            name, est = cfg[task]
            m = clone(est)
            rec = dict(base, task=task, model_name=name)
            if task == "regression":
                m.fit(X[tr], y[tr]); p = m.predict(X[T])
                rec.update(rmse=float(np.sqrt(mean_squared_error(y[T], p))), mae=float(mean_absolute_error(y[T], p)),
                           r2=float(r2_score(y[T], p)), bias=float((p - y[T]).mean()),
                           spearman=float(stats.spearmanr(p, y[T]).statistic), sd_y_test=float(y[T].std()))
            else:
                if len(np.unique(ycls[tr])) < 2:
                    continue
                m.fit(X[tr], ycls[tr]); p = m.predict_proba(X[T])[:, 1]
                has_two_classes = len(np.unique(ycls[T])) > 1
                rec.update(roc_auc=float(roc_auc_score(ycls[T], p)) if has_two_classes else float("nan"),
                           pr_auc=float(average_precision_score(ycls[T], p)) if has_two_classes else float("nan"))
            records.append(rec)
    print(group, ep, "random_seed", seed, round(time.time() - t0), "s", flush=True)
json.dump(records, open(out, "w"))
