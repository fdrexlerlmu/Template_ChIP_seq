#!/usr/bin/env bash
set -euo pipefail

export PATH="$HOME/.local/bin:$PATH"
PYTHON_BIN="python3"

# ==============================================================================
# Flexibles deepTools Metagene-Profiling (Mit strengem Promotor- & Fold-Change-Filter)
# + Replikat-Streuungsband (SD / SEM) für alle Plots
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

BIGWIG_DIR="${PROJECT_ROOT}/results/bigwig"
DIFFBIND_DIR="${PROJECT_ROOT}/results/diffbind"
GTF_FILE="${SCRIPT_DIR}/annotation/dm6_genes.gtf"
BED_FILE="${SCRIPT_DIR}/annotation/dm6_genes.bed"
FILTERED_BED="${SCRIPT_DIR}/annotation/dm6_stringent_responsive_genes.bed"
OUTPUT_DIR="${SCRIPT_DIR}/results_deeptools"

# Streuungsband-Modus: "SD" (Standardabweichung) oder "SEM" (Standardfehler)
BAND_MODE="${BAND_MODE:-SD}"

mkdir -p "${OUTPUT_DIR}/matrix"
mkdir -p "${OUTPUT_DIR}/plots/combined_states_stringent"
mkdir -p "${OUTPUT_DIR}/plots/combined_temporal_stringent"
mkdir -p "${OUTPUT_DIR}/merged_bigwigs"
mkdir -p "${SCRIPT_DIR}/annotation"

echo "========================================================"
echo " deepTools Pipeline gestartet (Strenger Zielgen-Modus)"
echo " Projekt-Root: ${PROJECT_ROOT}"
echo " Streuungsband: ±1 ${BAND_MODE}"
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

# 2. Replikate mergen (optional/zusätzlich für separate Auswertungen)
echo "Analysiere und merge Replikate..."
CONDITIONS=$(for bw in "${BIGWIG_DIR}"/*.bw; do
    [ -e "$bw" ] || continue
    basename "$bw" | sed -E 's/_rep[0-9]+\.bw$//'
done | sort -u)

for cond in $CONDITIONS; do
    files=("${BIGWIG_DIR}/${cond}"_rep*.bw)
    output_bw="${OUTPUT_DIR}/merged_bigwigs/${cond}.bw"
    if [ -f "$output_bw" ] || [ -L "$output_bw" ]; then
        continue
    elif [ "${#files[@]}" -eq 1 ] && [ -f "${files[0]}" ]; then
        ln -sf "${files[0]}" "$output_bw"
    elif [ "${#files[@]}" -gt 1 ]; then
        bigwigAverage -b "${files[@]}" -o "$output_bw" -p 8
    fi
done

# ==============================================================================
# 3. FINALER SWEET-SPOT-FILTER (Ziel: ~300-400 Gene)
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

print(" -> Durchsuche DiffBind-CSVs nach CrebA Primärkontrast...")

if csv_files:
    for csv_file in sorted(csv_files):
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

            sig_df = df_diff[(df_diff[fdr_col] < 0.001) & (df_diff[fold_col] >= 0.75)]

            for _, row in sig_df.iterrows():
                c = str(row[chrom_col])
                if not c.startswith('chr'):
                    c = 'chr' + c
                target_peaks.append((c, int(row[start_col]), int(row[end_col])))

            print(f"      -> {len(sig_df)} aktive Refed-Peaks extrahiert.")
        except Exception as e:
            print(f"Fehler beim Verarbeiten von {csv_file}: {e}")

print(f" -> Insgesamt {len(target_peaks)} Peaks gesammelt.")

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
            if g_chrom == p_chrom and max(g_start, p_start) < min(g_end, p_end):
                matched_lines.add(g_line)
                break

if matched_lines:
    print(f" -> EXAKT {len(matched_lines)} Gene entsprechen nun dem finalen Kern-Set.")
    final_lines = sorted(matched_lines)
else:
    print(" -> FEHLER: 0 Treffer nach Gen-Abgleich.")
    sys.exit(1)

with open(output_bed_path, "w") as out:
    out.writelines(final_lines)
EOF

# ==============================================================================
# 4. HILFSFUNKTION FÜR METAGENE-MATRIX & PLOT MIT REPLIKAT-STREUUNGSBAND
# ==============================================================================
run_metagene_with_band() {
    local matrix_file="$1"
    local output_pdf="$2"
    local plot_title="$3"
    shift 3
    local bw_files=("$@")

    if [ ${#bw_files[@]} -eq 0 ]; then
        return 0
    fi

    computeMatrix scale-regions \
      -S "${bw_files[@]}" \
      -R "${FILTERED_BED}" \
      --beforeRegionStartLength 2000 \
      --regionBodyLength 3000 \
      --afterRegionStartLength 2000 \
      --skipZeros \
      -out "${matrix_file}" \
      -p 8

    "$PYTHON_BIN" - "${matrix_file}" "${output_pdf}" "${plot_title}" "${BAND_MODE}" << 'EOF'
import os
import re
import sys
import gzip
import json
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

matrix_path, output_path, title, band_mode = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

with gzip.open(matrix_path, 'rt') as f:
    line = f.readline()
    header = json.loads(line[1:]) if line.startswith('@') else {}

sample_labels = header.get('sample_labels', [])
df = pd.read_csv(matrix_path, sep='\t', comment='@', header=None)

meta_cols = 6
total_data_cols = df.shape[1] - meta_cols
num_samples = len(sample_labels) if sample_labels else 1
bins_per_sample = total_data_cols // num_samples

cond_groups = defaultdict(list)
for i in range(num_samples):
    raw_label = sample_labels[i] if i < len(sample_labels) else f"Sample_{i+1}"
    base_name = os.path.basename(raw_label).replace(".bw", "")
    clean_name = re.sub(r'^(CrebA|Pol2)_ChIP_', '', base_name)
    cond_name = re.sub(r'_rep[0-9]+$', '', clean_name)

    start_col = meta_cols + i * bins_per_sample
    end_col = start_col + bins_per_sample
    rep_mat = df.iloc[:, start_col:end_col].apply(pd.to_numeric, errors='coerce').to_numpy(dtype=float)
    cond_groups[cond_name].append(rep_mat)

x = np.arange(bins_per_sample)
colors = ['#e41a1c', '#377eb8', '#4daf4a', '#ff7f00', '#984ea3', '#a65628', '#f781bf']

plt.figure(figsize=(9, 6))

for idx, (cond_name, rep_matrices) in enumerate(cond_groups.items()):
    stack = np.stack(rep_matrices, axis=0)  # (n_reps, n_genes, n_bins)
    n_reps = stack.shape[0]
    color = colors[idx % len(colors)]

    rep_mean_per_gene = np.nanmean(stack, axis=0)
    mean_profile = np.nanmean(rep_mean_per_gene, axis=0)

    if n_reps > 1:
        rep_sd_per_gene = np.nanstd(stack, axis=0, ddof=1)
        band_profile = np.nanmean(rep_sd_per_gene, axis=0)
        if band_mode.upper() == "SEM":
            band_profile = band_profile / np.sqrt(n_reps)
    else:
        band_profile = np.zeros_like(mean_profile)

    plt.plot(x, mean_profile, label=f"{cond_name} (n={n_reps})", color=color, linewidth=2)
    if n_reps > 1:
        plt.fill_between(
            x,
            mean_profile - band_profile,
            mean_profile + band_profile,
            color=color,
            alpha=0.22,
            linewidth=0
        )

tss_idx = int(bins_per_sample * (2000 / 7000))
tes_idx = int(bins_per_sample * (5000 / 7000))

plt.axvline(x=tss_idx, color='gray', linestyle='--', alpha=0.7, label='TSS')
plt.axvline(x=tes_idx, color='gray', linestyle=':', alpha=0.7, label='TES')
plt.xticks([0, tss_idx, tes_idx, bins_per_sample], ["-2.0 kb", "TSS", "TES", "+2.0 kb"])

plt.title(f"{title} (±1 {band_mode.upper()})", fontsize=14, fontweight='bold')
plt.xlabel("Genomic Region", fontsize=12)
plt.ylabel("Mean RPM Coverage", fontsize=12)
plt.legend(frameon=True, facecolor='white', edgecolor='none', fontsize=10)
plt.tight_layout()

plt.savefig(output_path, dpi=300)
plt.close()
print(f"Plot gespeichert unter: {output_path}")
EOF
}

# ==============================================================================
# 5. ALLE PLOTS ERZEUGEN (4h ZUSTÄNDE + ZEITVERLÄUFE)
# ==============================================================================
for target in Creba Pol2; do
    # 5a) 4h Zustände (alle Replikate)
    shopt -s nullglob
    state_files=("${BIGWIG_DIR}/${target}_ChIP_"*"_4h_rep"*.bw)
    shopt -u nullglob

    if [ ${#state_files[@]} -gt 0 ]; then
        echo "Berechne Matrix & Streuungsband-Plot für ${target} (4h Zustände)..."
        run_metagene_with_band \
            "${OUTPUT_DIR}/matrix/${target}_states_4h_reps_matrix.gz" \
            "${OUTPUT_DIR}/plots/combined_states_stringent/${target}_states_overlay_4h_band.pdf" \
            "${target} Target States (4h)" \
            "${state_files[@]}"
    else
        echo "Hinweis: Keine 4h-Replikat-BigWigs für ${target} gefunden."
    fi

    # 5b) Zeitlicher Verlauf über alle Zeitpunkte/Zustände (alle Replikate)
    shopt -s nullglob
    temporal_files=("${BIGWIG_DIR}/${target}_ChIP_"*"_rep"*.bw)
    shopt -u nullglob

    if [ ${#temporal_files[@]} -gt 0 ]; then
        echo "Berechne Matrix & Streuungsband-Plot für ${target} (Zeitverlauf / alle Replikate)..."
        run_metagene_with_band \
            "${OUTPUT_DIR}/matrix/${target}_temporal_all_reps_matrix.gz" \
            "${OUTPUT_DIR}/plots/combined_temporal_stringent/${target}_temporal_overlay_band.pdf" \
            "${target} Target Temporal Profile" \
            "${temporal_files[@]}"
    else
        echo "Hinweis: Keine zeitlichen Replikat-BigWigs für ${target} gefunden."
    fi
done

echo "========================================================"
echo " Fertig! Alle verfügbaren Replikat-Plots wurden erzeugt."
echo " Modus: ±1 ${BAND_MODE} (änderbar via BAND_MODE=SEM ./skript.sh)"
echo "========================================================"