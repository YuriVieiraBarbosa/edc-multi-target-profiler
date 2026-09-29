#!/usr/bin/env python
"""Executa predict_edc2.py (sem alteracoes), gera os graficos e guarda tudo em
uma pasta propria dentro de resultados/.

Aceita os mesmos argumentos de predict_edc2.py:
    python rodar_edc.py --input demo_edcs.csv --name-col nome --out demo

Cada execucao cria:  resultados/AAAAMMDD_HHMMSS_<out>/
    <out>_resumo.csv, <out>_detalhado.csv, <out>_rejeitadas.csv (se houver),
    <out>_dashboard.png, <out>_fracao_ativa.png,
    entrada.csv (copia do arquivo de entrada) e LEIAME.txt (comando, data, versoes).
"""
import argparse
import datetime
import os
import platform
import shutil
import subprocess
import sys

aqui = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(add_help=False)
ap.add_argument("--out", default="predicoes_consenso")
ap.add_argument("--input")
a, resto = ap.parse_known_args()

if "-h" in resto or "--help" in resto:
    sys.exit(subprocess.call([sys.executable, os.path.join(aqui, "predict_edc2.py")] + resto))

nome = os.path.basename(os.path.normpath(a.out)) or "predicoes_consenso"
agora = datetime.datetime.now()
base_res = os.path.join(os.getcwd(), "resultados")
os.makedirs(base_res, exist_ok=True)
pasta, k = None, 0
while pasta is None:
    cand = os.path.join(base_res, f"{agora:%Y%m%d_%H%M%S}_{nome}" + (f"_{k}" if k else ""))
    try:
        os.makedirs(cand)
        pasta = cand
    except FileExistsError:
        k += 1
prefixo = os.path.join(pasta, nome)

cmd = [sys.executable, os.path.join(aqui, "predict_edc2.py"), "--out", prefixo]
if a.input:
    cmd += ["--input", a.input]
cmd += resto
rc = subprocess.call(cmd)
if rc != 0:
    if not os.listdir(pasta):
        os.rmdir(pasta)
    sys.exit(rc)

if a.input and os.path.exists(a.input):
    shutil.copy2(a.input, os.path.join(pasta, "entrada.csv"))

sys.path.insert(0, aqui)
import plot_edc_results
plot_edc_results.gerar_todas(prefixo)


def _versao(mod):
    try:
        return __import__(mod).__version__
    except Exception:
        return "?"


linhas = [
    f"Execucao: {agora:%Y-%m-%d %H:%M:%S}",
    "Comando: python rodar_edc.py " + " ".join(sys.argv[1:]),
    f"Pasta de trabalho: {os.getcwd()}",
    f"Python {platform.python_version()} | rdkit {_versao('rdkit')} | "
    f"scikit-learn {_versao('sklearn')} | pandas {_versao('pandas')} | numpy {_versao('numpy')}",
    "",
    "Arquivos:",
] + [f"  {f}" for f in sorted(os.listdir(pasta))]
with open(os.path.join(pasta, "LEIAME.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(linhas) + "\n")

print(f"\nResultados desta execucao: {os.path.relpath(pasta)}")
for f in sorted(os.listdir(pasta)):
    print(f"   - {f}")
