#!/usr/bin/env python
"""Atalho: roda predict_edc2.py e, em seguida, a reagregação em dois eixos do consenso_v2.py.

Não altera predict_edc2.py nem nenhum modelo. Apenas chama os dois na ordem certa, com os cuidados
que o predict_edc2.py exige: ele identifica a camada pelo NOME da pasta de modelos, então precisa ser
executado a partir da raiz do pipeline e com os caminhos de modelo padrão. Passar `--model-dirs` com
caminho absoluto faz a camada Tox21 ser rotulada como ChEMBL e corrompe os votos.

Uso:
    python consenso_v2/rodar_consenso_v2.py --input moleculas.csv --name-col nome --out resultados/run1

Gera:
    <out>_detalhado.csv, <out>_resumo.csv          (do predict_edc2.py, formato original)
    <out>_v2_endpoints.csv, <out>_v2_resumo.csv,
    <out>_v2_cobertura.csv                          (do consenso_v2.py, dois eixos)
"""
import argparse
import os
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)


# Run predict_edc2.py and then the two-axis re-aggregation, in that order and from the pipeline root.
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="CSV de entrada com a coluna de SMILES")
    ap.add_argument("--smiles-col", default="smiles")
    ap.add_argument("--name-col", default=None)
    ap.add_argument("--out", default="predicoes_consenso", help="prefixo das saídas (pode ser caminho absoluto)")
    a = ap.parse_args()

    entrada = os.path.abspath(a.input)
    saida = os.path.abspath(a.out)
    os.makedirs(os.path.dirname(saida) or ".", exist_ok=True)

    cmd = [sys.executable, "predict_edc2.py", "--input", entrada, "--smiles-col", a.smiles_col, "--out", saida]
    if a.name_col:
        cmd += ["--name-col", a.name_col]
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    # cwd = raiz do pipeline: os caminhos models_final/ models_layers/ models_tox21/ são relativos e o
    # nome da pasta é o que define a camada de cada modelo.
    subprocess.run(cmd, cwd=RAIZ, env=env, check=True)

    det = f"{saida}_detalhado.csv"
    if not os.path.exists(det):
        sys.exit(f"predict_edc2.py não gerou {det}")
    subprocess.run([sys.executable, os.path.join(AQUI, "consenso_v2.py"), det, f"{saida}_v2"], check=True)


if __name__ == "__main__":
    main()
