"""Featurization: counted ECFP4 (2048 bits) + 2D RDKit descriptors."""
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
DESC_NAMES = [n for n, _ in Descriptors._descList if not n.startswith("fr_")]
_DESC_FUNCS = [f for n, f in Descriptors._descList if n in set(DESC_NAMES)]


# Parse a SMILES into an RDKit molecule, or None when parsing fails.
def mol_from_smiles(s):
    return Chem.MolFromSmiles(s)


# 2048-dimension ECFP4 count fingerprint, used as model input.
def ecfp4_counts(mol):
    fp = FPGEN.GetCountFingerprintAsNumPy(mol)
    return fp.astype(np.float32)


# 2048-bit binary ECFP4 fingerprint, used only for Tanimoto similarity and the applicability domain.
def ecfp4_bits(mol):
    return FPGEN.GetFingerprintAsNumPy(mol).astype(np.uint8)


# The 132 RDKit 2D descriptors; failures and non-finite values become zero.
def descriptors(mol):
    out = np.empty(len(_DESC_FUNCS), dtype=np.float64)
    for i, f in enumerate(_DESC_FUNCS):
        try:
            out[i] = f(mol)
        except Exception:
            out[i] = np.nan
    return out


# Build the 2180-column feature matrix (2048 fingerprint counts + 132 descriptors) and flag unparsable SMILES.
def featurize(smiles_list):
    mols = [mol_from_smiles(s) for s in smiles_list]
    ok = np.array([m is not None for m in mols])
    fp = np.vstack([ecfp4_counts(m) for m in np.array(mols, dtype=object)[ok]])
    de = np.vstack([descriptors(m) for m in np.array(mols, dtype=object)[ok]])
    de = np.nan_to_num(de, nan=0.0, posinf=0.0, neginf=0.0)
    bits = np.vstack([ecfp4_bits(m) for m in np.array(mols, dtype=object)[ok]])
    names = [f"ECFP4_{i}" for i in range(fp.shape[1])] + DESC_NAMES
    return np.hstack([fp, de]).astype(np.float32), bits, names, ok


# Bemis-Murcko scaffold as SMILES; returns an empty string for acyclic molecules, which therefore share one pseudo-group.
def murcko_scaffold(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)
    except Exception:
        return ""


# Split by Bemis-Murcko scaffold. NOTE: the code fills the TEST set with the largest groups, contrary to the docstring below.
def scaffold_split(smiles, test_size=0.2, seed=0, allow_random=False):
    """Split by Bemis-Murcko scaffold.

    ATTENTION: scaffolds of the LARGEST group go to the TEST set, filled greedily from
    largest to smallest up to `test_size`. The external test set is, therefore, extrapolation to the
    dominant chemical families: pessimistic and high-variance. This is the convention under which
    the deployed models were assessed, and so it was kept.

    Molecules without a scaffold (acyclic, or RDKit failure) each receive their own
    group. Grouping them under the empty scaffold would treat structurally unrelated compounds as
    a single chemical family, which already made the test degenerate into two Tox21 assays.

    If the greedy test does not reach half of the target, the function raises an error instead of silently
    falling back to random splitting: a random split labeled as "by scaffold"
    produces optimistic metrics under the wrong name. Pass allow_random=True to opt
    explicitly for that.
    """
    from collections import defaultdict
    groups = defaultdict(list)
    for i, s in enumerate(smiles):
        sc = murcko_scaffold(s)
        groups[sc if sc else ("__no_scaffold__", i)].append(i)
    sets = sorted(groups.values(), key=lambda x: (-len(x), x[0]))
    n_test = int(round(test_size * len(smiles)))
    train, test = [], []
    for grp in sets:
        if len(test) + len(grp) <= n_test:
            test.extend(grp)
        else:
            train.extend(grp)
    if len(test) < 0.5 * n_test:
        if not allow_random:
            raise ValueError(
                f"insufficient scaffold split: {len(test)} compounds in test versus "
                f"target of {n_test} ({len(sets)} scaffolds in {len(smiles)} compounds). "
                "Too few scaffolds for the requested fraction. Reduce test_size, group by another "
                "criterion, or pass allow_random=True to accept a RANDOM split "
                "(which is not a scaffold split and gives optimistic metrics).")
        rng = np.random.default_rng(seed)
        idx = rng.permutation(len(smiles))
        test = list(idx[:n_test]); train = list(idx[n_test:])
    return np.array(sorted(train)), np.array(sorted(test))


# Highest Tanimoto similarity of each query compound to the training set, computed in blocks.
def tanimoto_max_to_train(bits_query, bits_train, block=512):
    """Highest Tanimoto similarity of each query to the training set."""
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
