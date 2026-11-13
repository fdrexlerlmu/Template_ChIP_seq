
#!/usr/bin/env bash
set -uo pipefail

# Basis-Pfade (relativ zum Ordner Downstream_analysis/)
DATA_DIR="../results/diffbind"
BASE_ANNO_DIR="../results/annotation"
BASE_PLOT_DIR="../results/plots"
BASE_MOTIF_DIR="../results/motifs"
BASE_ENRICH_DIR="../results/enrichment"

echo "=== Starte Downstream-Analyse Pipeline (Multi-Target / Multi-Kontrast) ==="

# --------------------------------------------------
# DiffBind-Dateien rekursiv finden (Bash 3 kompatibel)
# --------------------------------------------------

csv_files=()

while IFS= read -r file; do
    csv_files[${#csv_files[@]}]="$file"
done < <(
    find "$DATA_DIR" -type f -name "differential_peaks_*.csv"
)

if [ ${#csv_files[@]} -eq 0 ]; then
    echo "[ERROR] Keine differential_peaks_*.csv Dateien in $DATA_DIR gefunden!"
    echo "Bitte zuerst Snakemake ausführen."
    exit 1
fi

echo "[INFO] ${#csv_files[@]} DiffBind-Dateien gefunden."

# ---------------------------------------------------------
# 2. Verarbeitung jedes Targets und Kontrasts
# ---------------------------------------------------------

for csv in "${csv_files[@]}"; do

    filename=$(basename "$csv")

    # Target aus übergeordnetem Ordner extrahieren
    target=$(basename "$(dirname "$csv")")

    # Kontrast aus Dateinamen extrahieren
    contrast="${filename#differential_peaks_}"
    contrast="${contrast%.csv}"

    # Target-spezifische UND kontrast-spezifische Ordner
    ANNO_DIR="${BASE_ANNO_DIR}/${target}/${contrast}"
    PLOT_DIR="${BASE_PLOT_DIR}/${target}/${contrast}"
    MOTIF_DIR="${BASE_MOTIF_DIR}/${target}/${contrast}"
    ENRICH_DIR="${BASE_ENRICH_DIR}/${target}/${contrast}"

    echo "--------------------------------------------------------"
    echo "[INFO] Target:   ${target}"
    echo "[INFO] Kontrast: ${contrast}"
    echo "[INFO] CSV:      ${csv}"
    echo "--------------------------------------------------------"

    # -----------------------------------------------------
    # Skip-Check
    # -----------------------------------------------------

    if [[ -f "${ANNO_DIR}/peakAnno.rds" ]] && \
       [[ -f "${MOTIF_DIR}/streme_up.rds" ]] && \
       [[ -f "${PLOT_DIR}/volcano_data.rds" ]]; then

        echo "[SKIP] ${target} / ${contrast} bereits analysiert."
        continue
    fi

    mkdir -p "$ANNO_DIR" "$PLOT_DIR" "$MOTIF_DIR" "$ENRICH_DIR"

    # -----------------------------------------------------
    # 1. Peak Annotation
    # -----------------------------------------------------

    echo "  -> [1/4] Peak Annotation..."

    Rscript scripts/01_annotation.R \
        "$csv" \
        "$ANNO_DIR"

    # -----------------------------------------------------
    # 2. Motif Enrichment
    # -----------------------------------------------------

    echo "  -> [2/4] Motif Enrichment (STREME)..."

    Rscript scripts/02_motifs.R \
        "$ANNO_DIR" \
        "$MOTIF_DIR"

    # -----------------------------------------------------
    # 3. Venn & Volcano
    # -----------------------------------------------------

    echo "  -> [3/4] Venn & Volcano..."

    Rscript scripts/03_venn_volcano.R \
        "$ANNO_DIR" \
        "$PLOT_DIR"

    # -----------------------------------------------------
    # 4. GO & KEGG Enrichment
    # -----------------------------------------------------

    echo "  -> [4/4] GO & KEGG Enrichment..."

    Rscript scripts/04_enrichment.R \
        "$PLOT_DIR" \
        "$ENRICH_DIR"

    echo "[SUCCESS] ${target} / ${contrast} abgeschlossen."

done

# ---------------------------------------------------------
# 3. Globale UpSet-Analyse
# ---------------------------------------------------------

echo "--------------------------------------------------------"
echo "[5/7] Bereite globale UpSet-Daten vor..."
echo "--------------------------------------------------------"

Rscript scripts/05_upset.R \
    "$BASE_PLOT_DIR" \
    "$BASE_ENRICH_DIR"

# ---------------------------------------------------------
# 4. Kinetische Profil-Clusteranalyse
# ---------------------------------------------------------

echo "--------------------------------------------------------"
echo "[6/7] Kinetische Profil-Clusteranalyse..."
echo "--------------------------------------------------------"

Rscript scripts/06_kinetic_clustering.R \
    "$BASE_PLOT_DIR" \
    "$BASE_ENRICH_DIR"

# ---------------------------------------------------------
# 5. Quarto-Bericht
# ---------------------------------------------------------

echo "--------------------------------------------------------"
echo "[7/7] Rendere globalen Quarto HTML-Bericht..."
echo "--------------------------------------------------------"

quarto render report_downstream.qmd

echo "========================================================"
echo "Alle Downstream-Analysen erfolgreich abgeschlossen!"
echo "========================================================"
