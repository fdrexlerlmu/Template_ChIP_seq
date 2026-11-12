#!/usr/bin/env bash
set -e

# Pfade relativ zum Skriptort (im Unterordner 'Downstream analysis'):
# '..' bringt dich hoch ins Hauptverzeichnis, wo 'results/' liegt.
DATA_DIR="../results/diffbind"
ANNO_DIR="../results/annotation"
PLOT_DIR="../results/plots"
MOTIF_DIR="../results/motifs"

INPUT_CSV="${DATA_DIR}/differential_peaks.csv"
ANNOTATED_CSV="${ANNO_DIR}/annotated_peaks.csv"

echo "=== Starte Downstream-Analyse Pipeline (aus 'Downstream analysis/') ==="

# 0. Verzeichnisse im Hauptordner erstellen (falls noch nicht vorhanden)
mkdir -p "${ANNO_DIR}" "${PLOT_DIR}" "${MOTIF_DIR}"

# 1. Peak Annotation
echo "[1/5] Annotiere Peaks..."
Rscript scripts/01_annotation.R "${INPUT_CSV}" "${ANNO_DIR}"

# 2. Motif Enrichment (STREME)
echo "[2/5] Führe Motif-Enrichment aus (STREME)..."
Rscript scripts/02_motifs.R "${ANNOTATED_CSV}" "${MOTIF_DIR}"

# 3. Venn- & Volcano-Datensätze
echo "[3/5] Bereite Venn & Volcano Datensätze vor..."
Rscript scripts/03_venn_volcano.R "${ANNOTATED_CSV}" "${PLOT_DIR}"

# 4. Funktionelle Anreicherungsanalyse (clusterProfiler: GO BP, GO CC, KEGG)
echo "[4/5] Führe GO & KEGG Enrichment-Analyse aus..."
Rscript scripts/04_enrichment.R "${ANNOTATED_CSV}" "${PLOT_DIR}"

# 5. Quarto-Bericht rendern
echo "[5/5] Rendere Quarto HTML-Bericht..."
quarto render report_downstream.qmd

echo "=== Downstream-Analyse erfolgreich abgeschlossen! ==="