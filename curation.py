"""Structural and bioactivity curation for the endocrine disruptor dataset.

RDKit standardization: cleanup -> parent fragment (salt/solvent removal) ->
charge neutralization -> tautomer canonicalization -> InChIKey.
"""
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

ORGANIC = {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Se", "Br", "I"}

_norm = rdMolStandardize.Normalizer()
_lfc = rdMolStandardize.LargestFragmentChooser()
_unch = rdMolStandardize.Uncharger()
_te = rdMolStandardize.TautomerEnumerator()
_te.SetMaxTautomers(200)


# Standardise one SMILES (salt stripped, neutralised, canonical tautomer) and return the InChIKey or the rejection reason.
def standardize(smiles):
    """Returns dict with standardized SMILES, InChIKey and the rejection reason."""
    res = {"std_smiles": None, "inchikey": None, "mw": np.nan,
           "n_frag_orig": np.nan, "undefined_stereo": np.nan, "reject_reason": None}
    if not isinstance(smiles, str) or not smiles.strip():
        res["reject_reason"] = "no_smiles"
        return res
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        res["reject_reason"] = "invalid_smiles"
        return res
    res["n_frag_orig"] = len(Chem.GetMolFrags(mol))
    try:
        mol = rdMolStandardize.Cleanup(mol)
        mol = _lfc.choose(mol)
        mol = _norm.normalize(mol)
        mol = _unch.uncharge(mol)
        mol = _te.Canonicalize(mol)
    except Exception as exc:
        res["reject_reason"] = f"standardization_failure:{type(exc).__name__}"
        return res
    if mol is None or mol.GetNumHeavyAtoms() == 0:
        res["reject_reason"] = "empty_after_standardization"
        return res
    els = {a.GetSymbol() for a in mol.GetAtoms()}
    if not els <= ORGANIC:
        res["reject_reason"] = "non_organic_elements:" + ",".join(sorted(els - ORGANIC))
        return res
    if mol.GetNumHeavyAtoms() < 6:
        res["reject_reason"] = "very_small"
        return res
    mw = Descriptors.MolWt(mol)
    if not (100.0 <= mw <= 1000.0):
        res["reject_reason"] = f"mw_out_of_range:{mw:.1f}"
        return res
    res["mw"] = mw
    try:
        res["std_smiles"] = Chem.MolToSmiles(mol)
        res["inchikey"] = Chem.MolToInchiKey(mol)
    except Exception as exc:
        res["std_smiles"] = None
        res["inchikey"] = None
        res["reject_reason"] = f"serialization_failure:{type(exc).__name__}"
        return res
    if not res["inchikey"]:
        res["reject_reason"] = "inchikey_empty"
        return res
    try:
        res["undefined_stereo"] = len(Chem.FindMolChiralCenters(
            mol, includeUnassigned=True, useLegacyImplementation=False)) - \
            len(Chem.FindMolChiralCenters(mol, includeUnassigned=False,
                                          useLegacyImplementation=False))
    except Exception:
        res["undefined_stereo"] = np.nan
    return res


# Standardise the unique SMILES in parallel; returns one row per distinct input structure.
def standardize_frame(smiles_series, n_jobs=8):
    from joblib import Parallel, delayed
    uniq = pd.Index(smiles_series.dropna().unique())
    recs = Parallel(n_jobs=n_jobs, batch_size=64)(delayed(standardize)(s) for s in uniq)
    out = pd.DataFrame(recs, index=uniq)
    out.index.name = "canonical_smiles"
    return out.reset_index()


# Convert a concentration in nM into its negative log10 (pActivity).
def to_pic50(value_nm):
    v = pd.to_numeric(value_nm, errors="coerce")
    v = v.where(v > 0)
    return 9.0 - np.log10(v)


# Median per group, discarding any group whose replicates span more than max_spread log units.
def aggregate_replicates(df, keys, value_col="pIC50", max_spread=1.0):
    """Median per group; discards groups with spread > max_spread log units."""
    g = df.groupby(keys, dropna=False)[value_col]
    agg = g.agg(n_records="size", pIC50_median="median", pIC50_min="min",
                pIC50_max="max", pIC50_std="std").reset_index()
    agg["spread"] = agg["pIC50_max"] - agg["pIC50_min"]
    agg["discordant"] = (agg["n_records"] > 1) & (agg["spread"] > max_spread)
    return agg
