"""Extração de atividades IC50 do dump SQLite do ChEMBL 37 para os alvos de
desregulação endócrina (AR, ERalfa, ERbeta, THRalfa, THRbeta, PPARgama).

Mantém proveniência completa por medida: composto, ensaio, alvo e referência.
"""
import sqlite3, sys, pandas as pd

ACC = {"AR": "P10275", "ESR1": "P03372", "ESR2": "Q92731",
       "PPARG": "P37231", "THRA": "P10827", "THRB": "P10828"}

TARGET_SQL = """
SELECT td.tid, td.chembl_id AS target_chembl_id, td.pref_name AS target_name,
       td.target_type, td.organism, cs.accession
FROM target_dictionary td
JOIN target_components tc ON td.tid = tc.tid
JOIN component_sequences cs ON tc.component_id = cs.component_id
WHERE cs.accession IN ({placeholders})
"""

ACT_SQL = """
SELECT
  act.activity_id, act.molregno,
  md.chembl_id              AS molecule_chembl_id,
  md.pref_name              AS molecule_name,
  md.max_phase, md.molecule_type, md.structure_type,
  cst.canonical_smiles, cst.standard_inchi_key,
  act.standard_type, act.standard_relation, act.standard_value, act.standard_units,
  act.pchembl_value, act.activity_comment, act.data_validity_comment,
  act.potential_duplicate, act.bao_endpoint,
  ass.chembl_id             AS assay_chembl_id,
  ass.description           AS assay_description,
  ass.assay_type, ass.assay_organism, ass.assay_tax_id, ass.confidence_score,
  ass.bao_format, ass.assay_test_type, ass.assay_category, ass.variant_id,
  ass.relationship_type, ass.cell_id, ass.tissue_id,
  td.chembl_id              AS target_chembl_id,
  td.pref_name              AS target_name,
  td.target_type, td.organism AS target_organism,
  d.chembl_id               AS document_chembl_id,
  d.doi, d.pubmed_id, d.year, d.journal, d.title AS document_title,
  src.src_short_name        AS activity_source,
  cp.mw_freebase, cp.alogp, cp.hba, cp.hbd, cp.psa, cp.rtb, cp.num_ro5_violations,
  cp.heavy_atoms, cp.aromatic_rings, cp.qed_weighted
FROM activities act
JOIN assays ass              ON act.assay_id = ass.assay_id
JOIN target_dictionary td    ON ass.tid = td.tid
JOIN molecule_dictionary md  ON act.molregno = md.molregno
LEFT JOIN compound_structures cst ON md.molregno = cst.molregno
LEFT JOIN compound_properties cp  ON md.molregno = cp.molregno
LEFT JOIN docs d             ON act.doc_id = d.doc_id
LEFT JOIN source src         ON act.src_id = src.src_id
WHERE td.tid IN ({tids})
  AND act.standard_type = 'IC50'
"""


def main(db_path, out_prefix="cache/chembl37"):
    con = sqlite3.connect(db_path)
    ph = ",".join("?" * len(ACC))
    tg = pd.read_sql(TARGET_SQL.format(placeholders=ph), con, params=list(ACC.values()))
    inv = {v: k for k, v in ACC.items()}
    tg["target_label"] = tg["accession"].map(inv)
    tg.to_csv(f"{out_prefix}_targets.csv", index=False)

    single = tg[(tg.target_type == "SINGLE PROTEIN") & (tg.organism == "Homo sapiens")]
    tids = sorted(set(single.tid.tolist()))
    acts = pd.read_sql(ACT_SQL.format(tids=",".join("?" * len(tids))), con, params=tids)
    acts = acts.merge(single[["target_chembl_id", "accession", "target_label"]].drop_duplicates(),
                      on="target_chembl_id", how="left")
    acts.to_parquet(f"{out_prefix}_ic50_raw.parquet", index=False)
    con.close()
    return tg, acts


if __name__ == "__main__":
    tg, acts = main(sys.argv[1])
    print("targets encontrados:", len(tg))
    print(tg[["target_chembl_id", "target_name", "target_type", "organism", "target_label"]]
          .to_string(index=False))
    print("\natividades IC50:", acts.shape)
    print(acts.groupby(["target_label", "assay_type"]).size().unstack(fill_value=0).to_string())
