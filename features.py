"""Featurização: ECFP4 contado (2048 bits) + descritores 2D do RDKit."""
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
def scaffold_split(smiles, test_size=0.2, seed=0, permitir_aleatorio=False):
    """Split por scaffold de Bemis-Murcko.

    ATENÇÃO: os scaffolds de MAIOR grupo vão para o TESTE, preenchido de forma gulosa do
    maior para o menor até `test_size`. O teste externo é, portanto, extrapolação para as
    famílias químicas dominantes: pessimista e de alta variância. Esta é a convenção com que
    os modelos implantados foram avaliados, e por isso foi mantida.

    Moléculas sem scaffold (acíclicas, ou falha do RDKit) recebem cada uma o seu próprio
    grupo. Agrupá-las sob o scaffold vazio trataria compostos sem relação estrutural como
    uma única família química, o que já tornou o teste degenerado em dois ensaios do Tox21.

    Se o teste guloso não alcançar metade do alvo, a função levanta erro em vez de cair em
    silêncio para divisão aleatória: uma divisão aleatória rotulada como "por scaffold"
    produz métricas otimistas sob um nome errado. Passe permitir_aleatorio=True para optar
    explicitamente por ela.
    """
    from collections import defaultdict
    groups = defaultdict(list)
    for i, s in enumerate(smiles):
        sc = murcko_scaffold(s)
        groups[sc if sc else ("__sem_scaffold__", i)].append(i)
    sets = sorted(groups.values(), key=lambda x: (-len(x), x[0]))
    n_test = int(round(test_size * len(smiles)))
    train, test = [], []
    for grp in sets:
        if len(test) + len(grp) <= n_test:
            test.extend(grp)
        else:
            train.extend(grp)
    if len(test) < 0.5 * n_test:
        if not permitir_aleatorio:
            raise ValueError(
                f"divisao por scaffold insuficiente: {len(test)} compostos no teste contra "
                f"alvo de {n_test} ({len(sets)} scaffolds em {len(smiles)} compostos). "
                "Poucos scaffolds para a fracao pedida. Reduza test_size, agrupe por outro "
                "criterio, ou passe permitir_aleatorio=True para aceitar divisao ALEATORIA "
                "(que nao e divisao por scaffold e da metricas otimistas).")
        rng = np.random.default_rng(seed)
        idx = rng.permutation(len(smiles))
        test = list(idx[:n_test]); train = list(idx[n_test:])
    return np.array(sorted(train)), np.array(sorted(test))


# Highest Tanimoto similarity of each query compound to the training set, computed in blocks.
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
