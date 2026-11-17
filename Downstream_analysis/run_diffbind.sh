#!/usr/bin/env bash
set -euo pipefail

mkdir -p logs/diffbind
echo "[INFO] Starte DESeq2 Analyse..."

# Ruft das saubere R-Skript auf und schreibt das Log
Rscript --vanilla scripts/run_diffbind.R 2>&1 | tee logs/diffbind/analysis.log

echo "[✓] Fertig!"