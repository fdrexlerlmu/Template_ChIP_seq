#!/usr/bin/env bash
set -euo pipefail

export PATH="$HOME/.local/bin:$PATH"
PYTHON_BIN="python3"

# ==============================================================================
# Flexibles deepTools Metagene-Profiling (Mit strengem Promotor- & Fold-Change-Filter)
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

BIGWIG_DIR="${PROJECT_ROOT}/results/bigwig"
DIFFBIND_DIR="${PROJECT_ROOT}/results/diffbind"
GTF_FILE="${SCRIPT_DIR}/annotation/dm6_genes.gtf"
BED_FILE="${SCRIPT_DIR}/annotation/dm6_genes.bed"
FILTERED_BED="${SCRIPT_DIR}/annotation/dm6_stringent_responsive_genes.bed"
OUTPUT_DIR="${SCRIPT_DIR}/results_deeptools"

mkdir -p "${OUTPUT_DIR}/matrix"
mkdir -p "${OUTPUT_DIR}/plots/combined_states_stringent"
mkdir -p "${OUTPUT_DIR}/plots/combined_temporal_stringent"
mkdir -p "${OUTPUT_DIR}/merged_bigwigs"
mkdir -p "${SCRIPT_DIR}/annotation"

echo "========================================================"
echo " deepTools Pipeline gestartet (Strenger Zielgen-Modus)"
echo " Projekt-Root: ${PROJECT_ROOT}"
echo "========================================================"

# 1. Basis-BED-Datei generieren
if [ ! -f "${BED_FILE}" ]; then
  echo "Generiere Basis-BED-Datei für kanonische Chromosomen..."
  "$PYTHON_BIN" - "${GTF_FILE}" "${BED_FILE}" << 'EOF'
import sys
gtf_path, bed_path = sys.argv[1], sys.argv[2]
VALID_CHROMS = {"chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"}
transcripts = {}
with open(gtf_path) as f:
    for line in f:
        if line.startswith("#") or not line.strip(): continue
        parts = line.strip().split("\t")
        if len(parts) < 9: continue
        raw_chrom = parts[0]
        chrom = raw_chrom if raw_chrom.startswith("chr") else "chr" + raw_chrom
        if chrom not in VALID_CHROMS: continue
        start, end, strand = int(parts[3]) - 1, int(parts[4]), parts[6]
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
EOF
fi

# 2. Replikate mergen
echo "Analysiere und merge Replikate..."
CONDITIONS=$(for bw in "${BIGWIG_DIR}"/*.bw; do
    basename "$bw" | sed -E 's/_rep[0-9]+\.bw$//'
done | sort -u)

for cond in $CONDITIONS; do
    files=("${BIGWIG_DIR}/${cond}"_rep*.bw)
    output_bw="${OUTPUT_DIR}/merged_bigwigs/${cond}.bw"
    if [ -f "$output_bw" ]; then
        continue
    elif [ "${#files[@]}" -eq 1 ] && [ -f "${files[0]}" ]; then
        ln -sf "${files[0]}" "$output_bw"
    elif [ "${#files[@]}" -gt 1 ]; then
        bigwigAverage -b "${files[@]}" -o "$output_bw" -p 8
    fi
done

# ==============================================================================
# 4. FINALER SWEET-SPOT-FILTER (Ziel: ~300-400 Gene)
# ==============================================================================
echo "Wende final abgestimmte Filter an (padj < 0.001, log2FC >= 0.75, TSS ± 100 bp)..."
"$PYTHON_BIN" - "${BED_FILE}" "${DIFFBIND_DIR}" "${FILTERED_BED}" << 'EOF'
import sys
import os
import glob
import pandas as pd

base_bed_path = sys.argv[1]
diffbind_dir = sys.argv[2]
output_bed_path = sys.argv[3]

target_peaks = []
csv_files = glob.glob(os.path.join(diffbind_dir, "**", "*.csv"), recursive=True)

print(f" -> Durchsuche DiffBind-CSVs nach CrebA Primärkontrast...")

if csv_files:
    for csv_file in csv_files:
        filename_lower = csv_file.lower()
        
        if 'creba' not in filename_lower:
            continue
        if 'refed_4h_vs_fasted_4h' not in filename_lower:
            continue
            
        print(f"   -> Lese Datei: {os.path.basename(csv_file)}")
        try:
            df_diff = pd.read_csv(csv_file)
            
            chrom_col = 'seqnames'
            start_col = 'start'
            end_col = 'end'
            fold_col = 'log2FoldChange'
            fdr_col = 'padj'
            
            df_diff[fdr_col] = pd.to_numeric(df_diff[fdr_col], errors='coerce')
            df_diff[fold_col] = pd.to_numeric(df_diff[fold_col], errors='coerce')
            df_diff = df_diff.dropna(subset=[fdr_col, fold_col])
            
            # Feintuning: padj < 0.001 und log2FC >= 0.75
            sig_df = df_diff[(df_diff[fdr_col] < 0.001) & (df_diff[fold_col] >= 0.75)]
            
            for _, row in sig_df.iterrows():
                c = str(row[chrom_col])
                if not c.startswith('chr'): c = 'chr' + c
                target_peaks.append((c, int(row[start_col]), int(row[end_col])))
                
            print(f"      -> {len(sig_df)} aktive Refed-Peaks extrahiert.")
        except Exception as e:
            print(f"Fehler beim Verarbeiten von {csv_file}: {e}")

print(f" -> Insgesamt {len(target_peaks)} Peaks gesammelt.")

# Basis-Gene einlesen und Promotor-Fenster auf TSS ± 100 bp setzen
gene_intervals = []
with open(base_bed_path) as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 6:
            chrom, start, end, strand = parts[0], int(parts[1]), int(parts[2]), parts[5]
            tss = start if strand == '+' else end
            gene_intervals.append((chrom, max(0, tss - 100), tss + 100, line))

matched_lines = set()
if target_peaks:
    for g_chrom, g_start, g_end, g_line in gene_intervals:
        for p_chrom, p_start, p_end in target_peaks:
            if g_chrom == p_chrom:
                if max(g_start, p_start) < min(g_end, p_end):
                    matched_lines.add(g_line)
                    break

if matched_lines:
    print(f" -> EXAKT {len(matched_lines)} Gene entsprechen nun dem finalen Kern-Set.")
    final_lines = list(matched_lines)
else:
    print(" -> FEHLER: 0 Treffer nach Gen-Abgleich.")
    sys.exit(1)

with open(output_bed_path, "w") as out:
    out.writelines(final_lines)
EOF

# ==============================================================================
# 4. PLOTS ERZEUGEN
# ==============================================================================
for target in Creba Pol2; do
    sub_files=(${OUTPUT_DIR}/merged_bigwigs/${target}_ChIP_*_4h.bw)
    
    if [ ${#sub_files[@]} -gt 0 ]; then
        matrix_file="${OUTPUT_DIR}/matrix/${target}_states_4h_stringent_matrix.gz"
        output_pdf="${OUTPUT_DIR}/plots/combined_states_stringent/${target}_states_overlay_4h_stringent.pdf"
        
        echo "Berechne stringente Matrix für ${target} (4h Zustände)..."
        computeMatrix scale-regions \
          -S "${sub_files[@]}" \
          -R "${FILTERED_BED}" \
          --beforeRegionStartLength 2000 \
          --regionBodyLength 3000 \
          --afterRegionStartLength 2000 \
          --skipZeros \
          -out "${matrix_file}" \
          -p 8
        
        echo "Generiere stringenten Overlay-Plot für ${target} (4h)..."
        "$PYTHON_BIN" - "${matrix_file}" "${output_pdf}" "${target} Target States (4h) " << 'EOF'
import sys
import gzip
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

matrix_path, output_path, title = sys.argv[1], sys.argv[2], sys.argv[3]

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
    label_name = raw_label.split("/")[-1].replace(".bw", "").replace("CrebA_ChIP_", "").replace("Pol2_ChIP_", "")
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
print(f"Stringenter Overlay-Plot gespeichert unter: {output_path}")
EOF
    fi
done

echo "========================================================"
echo " Fertig! Schau nach, wie viele Gene das Skript im Terminal"
echo " jetzt ausgibt (sollte deutlich näher an deiner Zielmenge sein)."
echo "========================================================"