"""Teste justo do efeito de pôr os grupos grandes de scaffold no treino.

Problema da comparação anterior: cada regra de split produz um TESTE diferente, então as métricas não são comparáveis.
Aqui o teste T é o mesmo nas três condições de treino:

  L = grupos grandes (os que features.scaffold_split manda para o teste: gulosos do maior ao menor, até 20%).
  S = o restante (grupos pequenos e scaffolds únicos).
  T = ~20% dos compostos, sorteados POR GRUPO de scaffold dentro de S (semente). R = S \ T.

  A  treino = R                      (grandes grupos fora do treino, como na composição de treino da regra atual)
  B  treino = R ∪ L                  (grandes grupos no treino, como na regra corrigida)
  C  treino = subamostra de R ∪ L com o tamanho de R (isola composição de quantidade de dados)

Algoritmo e hiperparâmetros: os escolhidos pelo protocolo no regime antigo (out/old_*.json), fixos nas três condições.
Não reexecuta a seleção por CV. Uso: cd work; python fair_test.py <ic50|layers|tox21> <endpoint> <semente> [<semente> ...]
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
grupo, ep = sys.argv[1], sys.argv[2]
seeds = [int(s) for s in sys.argv[3:]]
PASTA = {"ic50": ("data_ic50", "pIC50_final"), "layers": ("data_layers", "pAct_final"), "tox21": ("data_tox21", "ativo_especifico")}[grupo]
out = f"out/fair_{grupo}_{ep}.json"
if os.path.exists(out):
    sys.exit(0)

df = pd.read_csv(f"{PASTA[0]}/dataset_{ep}.csv")
cache = f"cache/{grupo}_{ep}.npz"
if os.path.exists(cache):
    z = np.load(cache, allow_pickle=True); X, bits, ok = z["X"], z["bits"], z["ok"]
else:
    X, bits, names, ok = F.featurize(df["std_smiles"].tolist())
    np.savez_compressed(cache, X=X, bits=bits, ok=ok)
df = df[ok].reset_index(drop=True)
sm = df["std_smiles"].tolist()
y = df[PASTA[1]].to_numpy(dtype=float)
ycls = y.astype(int) if grupo == "tox21" else (y >= 6.0).astype(int)
N = len(df)
sc = np.array([F.murcko_scaffold(s) for s in sm])

S, L = F.scaffold_split(sm, 0.2, 42)          # S = treino da regra atual; L = teste da regra atual (grupos grandes)
old = json.load(open(f"out/old_{grupo}_{ep}.json"))
res, rows = old["res"], old["rows"]


# Rebuild the algorithm and hyperparameters that the original protocol selected for this endpoint.
def config(tarefa):
    if grupo == "ic50":
        nome = res["regressao_melhor"] if tarefa == "regressao" else res.get("classificacao_melhor")
        chave = "melhores_params"
    elif grupo == "layers":
        nome = res["melhor_modelo"] if tarefa == "regressao" else res.get("melhor_classificador")
        chave = "params"
    else:
        nome, chave = res["melhor_modelo"], "params"
    if not nome:
        return None
    r = [x for x in rows if x["modelo"] == nome and (grupo == "tox21" or x["tarefa"] == tarefa)][0]
    zoo = (M.reg_zoo(X.shape[1]) if tarefa == "regressao" else (MT.zoo() if grupo == "tox21" else M.clf_zoo()))
    est = clone(zoo[nome][0]).set_params(**json.loads(r[chave]))
    if "n_jobs" in est.get_params():
        est.set_params(n_jobs=NJ)
    return nome, est


tarefas = ["classificacao"] if grupo == "tox21" else ["regressao", "classificacao"]
cfg = {t: config(t) for t in tarefas}
registros = []
t0 = time.time()
for seed in seeds:
    rng = np.random.default_rng(seed)
    grp = defaultdict(list)
    for i in S:
        grp[sc[i]].append(i)
    ch = list(grp.values())
    alvo = int(round(0.2 * N)); T = []
    for k in rng.permutation(len(ch)):
        if len(T) + len(ch[k]) <= alvo:
            T.extend(ch[k])
        if len(T) >= alvo:
            break
    T = np.array(sorted(T)); Tset = set(T.tolist())
    R = np.array([i for i in S if i not in Tset]); RL = np.concatenate([R, L])
    C = np.sort(np.random.default_rng(seed + 1000).choice(RL, size=len(R), replace=False))
    conds = {"A_sem_grandes": R, "B_com_grandes": RL, "C_tamanho_igual": C}
    for cn, tr in conds.items():
        simT = F.tanimoto_max_to_train(bits[T], bits[tr])
        base = dict(grupo=grupo, endpoint=ep, semente=seed, condicao=cn, n_treino=int(len(tr)), n_teste=int(len(T)),
                    frac_treino_de_L=float(np.isin(tr, L).mean()), scaffolds_teste=int(len(set(sc[T]))),
                    tanimoto_max_teste_mediana=float(np.median(simT)), prevalencia_teste=float(ycls[T].mean()),
                    prevalencia_treino=float(ycls[tr].mean()))
        for tarefa in tarefas:
            if cfg[tarefa] is None:
                continue
            nome, est = cfg[tarefa]
            m = clone(est)
            rec = dict(base, tarefa=tarefa, modelo=nome)
            if tarefa == "regressao":
                m.fit(X[tr], y[tr]); p = m.predict(X[T])
                rec.update(rmse=float(np.sqrt(mean_squared_error(y[T], p))), mae=float(mean_absolute_error(y[T], p)),
                           r2=float(r2_score(y[T], p)), vies=float((p - y[T]).mean()),
                           spearman=float(stats.spearmanr(p, y[T]).statistic), sd_y_teste=float(y[T].std()))
            else:
                if len(np.unique(ycls[tr])) < 2:
                    continue
                m.fit(X[tr], ycls[tr]); p = m.predict_proba(X[T])[:, 1]
                dois = len(np.unique(ycls[T])) > 1
                rec.update(roc_auc=float(roc_auc_score(ycls[T], p)) if dois else float("nan"),
                           pr_auc=float(average_precision_score(ycls[T], p)) if dois else float("nan"))
            registros.append(rec)
    print(grupo, ep, "semente", seed, round(time.time() - t0), "s", flush=True)
json.dump(registros, open(out, "w"))
