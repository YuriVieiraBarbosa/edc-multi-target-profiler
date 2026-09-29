"""Featurização: ECFP4 contado (2048 bits) + descritores 2D do RDKit."""
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
DESC_NAMES = [n for n, _ in Descriptors._descList if not n.startswith("fr_")]
_DESC_FUNCS = [f for n, f in Descriptors._descList if n in set(DESC_NAMES)]


def mol_from_smiles(s):
    return Chem.MolFromSmiles(s)


def ecfp4_counts(mol):
    fp = FPGEN.GetCountFingerprintAsNumPy(mol)
    return fp.astype(np.float32)


def ecfp4_bits(mol):
    return FPGEN.GetFingerprintAsNumPy(mol).astype(np.uint8)


def descriptors(mol):
    out = np.empty(len(_DESC_FUNCS), dtype=np.float64)
    for i, f in enumerate(_DESC_FUNCS):
        try:
            out[i] = f(mol)
        except Exception:
            out[i] = np.nan
    return out


def featurize(smiles_list):
    mols = [mol_from_smiles(s) for s in smiles_list]
    ok = np.array([m is not None for m in mols])
    fp = np.vstack([ecfp4_counts(m) for m in np.array(mols, dtype=object)[ok]])
    de = np.vstack([descriptors(m) for m in np.array(mols, dtype=object)[ok]])
    de = np.nan_to_num(de, nan=0.0, posinf=0.0, neginf=0.0)
    bits = np.vstack([ecfp4_bits(m) for m in np.array(mols, dtype=object)[ok]])
    names = [f"ECFP4_{i}" for i in range(fp.shape[1])] + DESC_NAMES
    return np.hstack([fp, de]).astype(np.float32), bits, names, ok


def murcko_scaffold(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)
    except Exception:
        return ""


def scaffold_split(smiles, test_size=0.2, seed=0):
    """Split por scaffold de Bemis-Murcko: scaffolds maiores vão para treino."""
    from collections import defaultdict
    groups = defaultdict(list)
    for i, s in enumerate(smiles):
        groups[murcko_scaffold(s)].append(i)
    sets = sorted(groups.values(), key=lambda x: (-len(x), x[0]))
    n_test = int(round(test_size * len(smiles)))
    train, test = [], []
    rng = np.random.default_rng(seed)
    for grp in sets:
        if len(test) + len(grp) <= n_test:
            test.extend(grp)
        else:
            train.extend(grp)
    if len(test) < 0.5 * n_test:  # fallback para datasets com poucos scaffolds
        idx = rng.permutation(len(smiles))
        test = list(idx[:n_test]); train = list(idx[n_test:])
    return np.array(sorted(train)), np.array(sorted(test))


def tanimoto_max_to_train(bits_query, bits_train, block=512):
    """Maior similaridade de Tanimoto de cada consulta ao conjunto de treino."""
    q = bits_query.astype(np.float32)
    t = bits_train.astype(np.float32)
    tn = t.sum(1)
    out = np.empty(len(q), dtype=np.float32)
    for i in range(0, len(q), block):
        qb = q[i:i + block]
        inter = qb @ t.T
        union = qb.sum(1)[:, None] + tn[None, :] - inter
        out[i:i + block] = np.max(inter / np.maximum(union, 1e-9), axis=1)
    return out
