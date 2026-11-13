#!/usr/bin/env bash
set -euo pipefail

# Basis-Pfade (relativ zum Ordner 'Downstream_analysis/')
DATA_DIR="../results/diffbind"
BASE_ANNO_DIR="../results/annotation"
BASE_PLOT_DIR="../results/plots"
BASE_MOTIF_DIR="../results/motifs"
BASE_ENRICH_DIR="../results/enrichment"

echo "=== Starte Downstream-Analyse Pipeline (Multi-Kontrast Modus) ==="

# Alle generierten DiffBind CSVs automatisch finden
shopt -s nullglob
csv_files=("${DATA_DIR}"/differential_peaks_*.csv)

if [ ${#csv_files[@]} -eq 0 ]; then
    echo "[ERROR] Keine differential_peaks_*.csv Dateien in ${DATA_DIR} gefunden! Bitte zuerst Snakemake ausführen."
    exit 1
fi

# Schleife über jeden gefundenen Kontrast
for csv in "${csv_files[@]}"; do
    # Kontrast-Namen aus dem Dateinamen extrahieren (z.B. "Fasted_4h_vs_Fed_4h")
    filename=$(basename "$csv")
    contrast="${filename#differential_peaks_}"
    contrast="${contrast%.csv}"

    echo "--------------------------------------------------------"
    echo "Verarbeite Kontrast: ${contrast}"
    echo "--------------------------------------------------------"

    # Kontrast-spezifische Unterordner definieren und anlegen
    ANNO_DIR="${BASE_ANNO_DIR}/${contrast}"
    PLOT_DIR="${BASE_PLOT_DIR}/${contrast}"
    MOTIF_DIR="${BASE_MOTIF_DIR}/${contrast}"
    ENRICH_DIR="${BASE_ENRICH_DIR}/${contrast}"

    mkdir -p "${ANNO_DIR}" "${PLOT_DIR}" "${MOTIF_DIR}" "${ENRICH_DIR}"

    # 1. Peak Annotation (Übergabe der konkreten CSV-Datei als Argument 1)
    echo "[1/5] Annotiere Peaks..."
    Rscript scripts/01_annotation.R "${csv}" "${ANNO_DIR}"

    # 2. Motif Enrichment (STREME)
    echo "[2/5] Führe Motif-Enrichment aus (STREME)..."
    Rscript scripts/02_motifs.R "${ANNO_DIR}" "${MOTIF_DIR}"

    # 3. Venn- & Volcano-Datensätze vorbereiten
    echo "[3/5] Bereite Venn & Volcano Datensätze vor..."
    Rscript scripts/03_venn_volcano.R "${ANNO_DIR}" "${PLOT_DIR}"

    # 4. Funktionelle Anreicherungsanalyse (GO & KEGG)
    echo "[4/5] Führe GO & KEGG Enrichment-Analyse aus..."
    Rscript scripts/04_enrichment.R "${PLOT_DIR}" "${ENRICH_DIR}"

    # 5. Quarto-Bericht für diesen spezifischen Kontrast rendern
    echo "[5/5] Rendere Quarto HTML-Bericht für ${contrast}..."
    quarto render report_downstream.qmd -P contrast:"${contrast}"
    
    # Bericht kontrast-spezifisch umbenennen, damit er nicht überschrieben wird
    if [ -f "report_downstream.html" ]; then
        mv report_downstream.html "report_downstream_${contrast}.html"
    fi
done

echo "=== Alle Downstream-Analysen erfolgreich abgeschlossen! ==="