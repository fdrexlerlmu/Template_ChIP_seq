#!/usr/bin/env bash
set -euo pipefail

# Ordner und Dateinamen definieren
RAW_DIR="rawdata"
OUTPUT_FILE="Metadata.txt"

if [ ! -d "$RAW_DIR" ]; then
    echo "[ERROR] Ordner '$RAW_DIR' nicht gefunden! Bitte im Projektverzeichnis anlegen." >&2
    exit 1
fi

# Temporäre Datei für die Proben-Erfassung
tmp_list=$(mktemp)

# Alle gzipped FASTQ/TXT-Dateien einlesen
for f in "$RAW_DIR"/*.{fastq,txt}.gz; do
    [ -e "$f" ] || continue
    filename=$(basename "$f")
    
    # Suffixe wie _R1, _R2, _1, _2 oder Dateiendungen sauber abschneiden, um die reine Sample_ID zu erhalten
    if [[ "$filename" =~ ^(.+)(_R[12]|_[12])\.(fastq|txt)\.gz$ ]]; then
        sample_id="${BASH_REMATCH[1]}"
        echo "$sample_id" >> "$tmp_list"
    elif [[ "$filename" =~ ^(.+)\.(fastq|txt)\.gz$ ]]; then
        sample_id="${BASH_REMATCH[1]}"
        echo "$sample_id" >> "$tmp_list"
    fi
done

# Vollständigen Header mit allen gewünschten Spalten schreiben
echo -e "Sample_ID\tTarget\tSample_Type\tInput_Sample_ID\tGenotype\tCondition\tBio_Replicate\tTissue\tAntibody\tFragment_Length\tSeq_Mode" > "$OUTPUT_FILE"

# Eindeutige Proben-IDs einlesen, sortieren und alle restlichen Spalten komplett leer lassen
while read -r sample_id; do
    # Jede Spalte ist mit einem Tabulator getrennt, aber leer – so kannst du alles händisch eintragen
    echo -e "${sample_id}\t\t\t\t\t\t\t\t\t\t"
done < <(sort -u "$tmp_list") >> "$OUTPUT_FILE"

# Aufräumen
rm -f "$tmp_list"

# Anzahl der geschriebenen Proben ermitteln (ohne Header)
num_samples=$(($(wc -l < "$OUTPUT_FILE") - 1))
echo "[SUCCESS] '$OUTPUT_FILE' wurde mit ${num_samples} Proben generiert!"