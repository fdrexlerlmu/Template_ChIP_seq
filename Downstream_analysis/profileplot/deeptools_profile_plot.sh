#!/usr/bin/env bash
set -euo pipefail

# PATH explizit erweitern, damit pipx-Befehle immer gefunden werden
export PATH="$HOME/.local/bin:$PATH"

# ==============================================================================
# Flexibles deepTools Metagene-Profiling (Perfekt skalierte Overlay-Plots)
# ==============================================================================

# 1. Dynamische Pfad-Ermittlung
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

BIGWIG_DIR="${PROJECT_ROOT}/results/bigwig"
GTF_FILE="${SCRIPT_DIR}/annotation/dm6_genes.gtf"
BED_FILE="${SCRIPT_DIR}/annotation/dm6_genes.bed"
OUTPUT_DIR="${SCRIPT_DIR}/results_deeptools"

mkdir -p "${OUTPUT_DIR}/matrix"
mkdir -p "${OUTPUT_DIR}/plots/combined_states"
mkdir -p "${OUTPUT_DIR}/plots/combined_temporal"
mkdir -p "${OUTPUT_DIR}/merged_bigwigs"
mkdir -p "${SCRIPT_DIR}/annotation"

echo "========================================================"
echo " deepTools Pipeline gestartet (Skalierungs-Fix)"
echo " Projekt-Root: ${PROJECT_ROOT}"
echo "========================================================"

# 2. Automatische Konvertierung von GTF zu BED (kanonische Chromosomen mit 'chr')
if [ ! -f "${BED_FILE}" ]; then
  echo "Generiere BED-Datei für kanonische Chromosomen..."
  python3 - "${GTF_FILE}" "${BED_FILE}" << 'EOF'
import sys

gtf_path = sys.argv[1]
bed_path = sys.argv[2]
VALID_CHROMS = {"chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"}

transcripts = {}
with open(gtf_path) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        
        raw_chrom = parts[0]
        chrom = raw_chrom if raw_chrom.startswith("chr") else "chr" + raw_chrom
        if chrom not in VALID_CHROMS:
            continue
            
        start = int(parts[3]) - 1
        end = int(parts[4])
        strand = parts[6]
        
        attr = parts[8]
        if "transcript_id" in attr:
            tid_start = attr.find("transcript_id") + 15
            tid_end = attr.find(";", tid_start)
            if tid_end == -1: tid_end = len(attr)
            tid = attr[tid_start:tid_end].strip().strip('"')
            
            if tid not in transcripts:
                transcripts[tid] = {"chrom": chrom, "start": start, "end": end, "strand": strand}
            else:
                transcripts[tid]["start"] = min(transcripts[tid]["start"], start)
                transcripts[tid]["end"] = max(transcripts[tid]["end"], end)

with open(bed_path, "w") as out:
    for tid, data in transcripts.items():
        out.write(f"{data['chrom']}\t{data['start']}\t{data['end']}\t{tid}\t0\t{data['strand']}\n")

print(f"Erfolgreich {len(transcripts)} Transkripte exportiert!")
EOF
fi

# 3. Replikate dynamisch gruppieren, cachen und parallel gemerged
echo "Analysiere und merge Replikate pro Zustand (Turbo-Modus)..."
CONDITIONS=$(for bw in "${BIGWIG_DIR}"/*.bw; do
    basename "$bw" | sed -E 's/_rep[0-9]+\.bw$//'
done | sort -u)

for cond in $CONDITIONS; do
    files=("${BIGWIG_DIR}/${cond}"_rep*.bw)
    num_files=${#files[@]}
    output_bw="${OUTPUT_DIR}/merged_bigwigs/${cond}.bw"
    
    if [ -f "$output_bw" ]; then
        echo " -> [Cache] Verwende bestehendes gemergtes BigWig für: ${cond}"
    elif [ "$num_files" -eq 1 ] && [ -f "${files[0]}" ]; then
        ln -sf "${files[0]}" "$output_bw"
    elif [ "$num_files" -gt 1 ]; then
        echo " -> [Berechnung] Mittle ${num_files} Replikate für ${cond} mit 8 Kernen..."
        bigwigAverage -b "${files[@]}" -o "$output_bw" -p 8
    fi
done

# ==============================================================================
# 4. ECHTE OVERLAY-PLOTS MIT EXakter PROPORTIONALER SKALIERUNG ERZEUGEN
# ==============================================================================

# A) ZUSTANDSVERGLEICH BEI 4H (Fed vs Fasted vs Refed in EINEM Plot)
for target in Creba Pol2; do
    sub_files=(${OUTPUT_DIR}/merged_bigwigs/${target}_ChIP_*_4h.bw)
    
    if [ ${#sub_files[@]} -gt 0 ]; then
        matrix_file="${OUTPUT_DIR}/matrix/${target}_states_4h_matrix.gz"
        output_pdf="${OUTPUT_DIR}/plots/combined_states/${target}_states_overlay_4h.pdf"
        
        echo "Berechne Matrix für ${target} (4h Zustände)..."
        computeMatrix scale-regions \
          -S "${sub_files[@]}" \
          -R "${BED_FILE}" \
          --beforeRegionStartLength 2000 \
          --regionBodyLength 3000 \
          --afterRegionStartLength 2000 \
          --skipZeros \
          -out "${matrix_file}" \
          -p 8
        
        echo "Generiere perfekt skalierten Overlay-Plot für ${target} (4h)..."
        python3 - "${matrix_file}" "${output_pdf}" "${target} Nutritional States (4h)" << 'EOF'
import sys
import gzip
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

matrix_path = sys.argv[1]
output_path = sys.argv[2]
title = sys.argv[3]

with gzip.open(matrix_path, 'rt') as f:
    line = f.readline()
    header = json.loads(line[1:]) if line.startswith('@') else {}

sample_labels = header.get('sample_labels', [])
df = pd.read_csv(matrix_path, sep='\t', comment='@', header=None)

meta_cols = 6
total_data_cols = df.shape[1] - meta_cols
num_samples = len(sample_labels) if sample_labels else 3
bins_per_sample = total_data_cols // num_samples

x = np.arange(bins_per_sample)
colors = ['#e41a1c', '#377eb8', '#4daf4a', '#ff7f00', '#984ea3']

plt.figure(figsize=(9, 6))

for i in range(num_samples):
    start_col = meta_cols + i * bins_per_sample
    end_col = start_col + bins_per_sample
    
    sample_data = df.iloc[:, start_col:end_col].apply(pd.to_numeric, errors='coerce')
    mean_profile = sample_data.mean(axis=0).values
    
    raw_label = sample_labels[i] if i < len(sample_labels) else f"Sample {i+1}"
    label_name = raw_label.split("/")[-1].replace(".bw", "").replace("Creba_ChIP_", "").replace("Pol2_ChIP_", "")
    
    plt.plot(x, mean_profile, label=label_name, color=colors[i % len(colors)], linewidth=2)

# Exakte proportionale Berechnung basierend auf den computeMatrix Parametern:
# Vor Start: 2000bp | Genkörper: 3000bp | Nach Ende: 2000bp (Gesamt: 7000bp)
tss_idx = int(bins_per_sample * (2000 / 7000))
tes_idx = int(bins_per_sample * (5000 / 7000))

plt.axvline(x=tss_idx, color='gray', linestyle='--', alpha=0.7, label='TSS')
plt.axvline(x=tes_idx, color='gray', linestyle=':', alpha=0.7, label='TES')
plt.xticks([0, tss_idx, tes_idx, bins_per_sample], ["-2.0 kb", "TSS", "TES", "+2.0 kb"])

plt.title(title, fontsize=14, fontweight='bold')
plt.xlabel("Genomic Region", fontsize=12)
plt.ylabel("Mean RPM Coverage", fontsize=12)
plt.legend(frameon=True, facecolor='white', edgecolor='none', fontsize=11)
plt.tight_layout()

plt.savefig(output_path, dpi=300)
plt.close()
print(f"Perfekter Overlay-Plot gespeichert unter: {output_path}")
EOF
    fi
done

# B) ZEITLICHER VERLAUF IM REFED-ZUSTAND (Refed 2h vs Refed 4h in EINEM Plot)
for target in Creba Pol2; do
    sub_files=(${OUTPUT_DIR}/merged_bigwigs/${target}_ChIP_Refed_*.bw)
    
    if [ ${#sub_files[@]} -gt 1 ]; then
        matrix_file="${OUTPUT_DIR}/matrix/${target}_temporal_refed_matrix.gz"
        output_pdf="${OUTPUT_DIR}/plots/combined_temporal/${target}_temporal_refed_overlay.pdf"
        
        echo "Berechne zeitliche Matrix für ${target} (Refed Verlauf)..."
        computeMatrix scale-regions \
          -S "${sub_files[@]}" \
          -R "${BED_FILE}" \
          --beforeRegionStartLength 2000 \
          --regionBodyLength 3000 \
          --afterRegionStartLength 2000 \
          --skipZeros \
          -out "${matrix_file}" \
          -p 8
        
        echo "Generiere zeitlichen Overlay-Plot für ${target} (Refed)..."
        python3 - "${matrix_file}" "${output_pdf}" "${target} Temporal Dynamics (Refed State)" << 'EOF'
import sys
import gzip
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

matrix_path = sys.argv[1]
output_path = sys.argv[2]
title = sys.argv[3]

with gzip.open(matrix_path, 'rt') as f:
    line = f.readline()
    header = json.loads(line[1:]) if line.startswith('@') else {}

sample_labels = header.get('sample_labels', [])
df = pd.read_csv(matrix_path, sep='\t', comment='@', header=None)

meta_cols = 6
total_data_cols = df.shape[1] - meta_cols
num_samples = len(sample_labels) if sample_labels else 2
bins_per_sample = total_data_cols // num_samples

x = np.arange(bins_per_sample)
colors = ['#377eb8', '#ff7f00', '#4daf4a']

plt.figure(figsize=(9, 6))

for i in range(num_samples):
    start_col = meta_cols + i * bins_per_sample
    end_col = start_col + bins_per_sample
    
    sample_data = df.iloc[:, start_col:end_col].apply(pd.to_numeric, errors='coerce')
    mean_profile = sample_data.mean(axis=0).values
    
    raw_label = sample_labels[i] if i < len(sample_labels) else f"Sample {i+1}"
    label_name = raw_label.split("/")[-1].replace(".bw", "").replace("Creba_ChIP_", "").replace("Pol2_ChIP_", "")
    
    plt.plot(x, mean_profile, label=label_name, color=colors[i % len(colors)], linewidth=2)

tss_idx = int(bins_per_sample * (2000 / 7000))
tes_idx = int(bins_per_sample * (5000 / 7000))

plt.axvline(x=tss_idx, color='gray', linestyle='--', alpha=0.7, label='TSS')
plt.axvline(x=tes_idx, color='gray', linestyle=':', alpha=0.7, label='TES')
plt.xticks([0, tss_idx, tes_idx, bins_per_sample], ["-2.0 kb", "TSS", "TES", "+2.0 kb"])

plt.title(title, fontsize=14, fontweight='bold')
plt.xlabel("Genomic Region", fontsize=12)
plt.ylabel("Mean RPM Coverage", fontsize=12)
plt.legend(frameon=True, facecolor='white', edgecolor='none', fontsize=11)
plt.tight_layout()

plt.savefig(output_path, dpi=300)
plt.close()
print(f"Zeitlicher Overlay-Plot gespeichert unter: {output_path}")
EOF
    fi
done

echo "========================================================"
echo " Fertig! Die perfekt skalierten Plots liegen unter:"
echo " ${OUTPUT_DIR}/plots/combined_states/"
echo " ${OUTPUT_DIR}/plots/combined_temporal/"
echo "========================================================"