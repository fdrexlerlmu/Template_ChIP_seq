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
mkdir -p "${OUTPUT_DIR}/merged_bigwigs"
mkdir -p "${SCRIPT_DIR}/annotation"

# Alte Plots ohne Streuungsband löschen
rm -f "${OUTPUT_DIR}/plots/combined_states_stringent/"*_stringent.pdf

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
        chrom = parts[0] if parts[0].startswith("chr") else "chr" + parts[0]
        if chrom not in VALID_CHROMS: continue
        start, end, strand = int(parts[3]) - 1, int(parts[4]), parts[6]
        attr = parts[8]
        if "transcript_id" in attr:
            i = attr.find("transcript_id") + 15
            j = attr.find(";", i)
            tid = attr[i:(j if j != -1 else len(attr))].strip().strip('"')
            if tid not in transcripts:
                transcripts[tid] = {"chrom": chrom, "start": start, "end": end, "strand": strand}
            else:
                transcripts[tid]["start"] = min(transcripts[tid]["start"], start)
                transcripts[tid]["end"] = max(transcripts[tid]["end"], end)
with open(bed_path, "w") as out:
    for tid, d in transcripts.items():
        out.write(f"{d['chrom']}\t{d['start']}\t{d['end']}\t{tid}\t0\t{d['strand']}\n")
EOF
fi

# 2. Replikate mergen
echo "Analysiere und merge Replikate..."
for cond in $(for bw in "${BIGWIG_DIR}"/*.bw; do basename "$bw" | sed -E 's/_rep[0-9]+\.bw$//'; done | sort -u); do
    files=("${BIGWIG_DIR}/${cond}"_rep*.bw)
    out_bw="${OUTPUT_DIR}/merged_bigwigs/${cond}.bw"
    [ -e "$out_bw" ] && continue
    if [ "${#files[@]}" -eq 1 ]; then
        ln -sf "${files[0]}" "$out_bw"
    elif [ "${#files[@]}" -gt 1 ]; then
        bigwigAverage -b "${files[@]}" -o "$out_bw" -p 8
    fi
done

# ==============================================================================
# 4. FINALER SWEET-SPOT-FILTER (Ziel: ~300-400 Gene)
# ==============================================================================
echo "Wende final abgestimmte Filter an (padj < 0.001, log2FC >= 0.75, TSS ± 100 bp)..."
"$PYTHON_BIN" - "${BED_FILE}" "${DIFFBIND_DIR}" "${FILTERED_BED}" << 'EOF'
import sys, os, glob
import pandas as pd

base_bed_path, diffbind_dir, output_bed_path = sys.argv[1], sys.argv[2], sys.argv[3]
target_peaks = []

for csv_file in sorted(glob.glob(os.path.join(diffbind_dir, "**", "*.csv"), recursive=True)):
    name = csv_file.lower()
    if 'creba' not in name or 'refed_4h_vs_fasted_4h' not in name:
        continue
    print(f"   -> Lese Datei: {os.path.basename(csv_file)}")
    df = pd.read_csv(csv_file)
    df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
    df['log2FoldChange'] = pd.to_numeric(df['log2FoldChange'], errors='coerce')
    sig = df.dropna(subset=['padj', 'log2FoldChange'])
    sig = sig[(sig['padj'] < 0.001) & (sig['log2FoldChange'] >= 0.75)]
    for _, r in sig.iterrows():
        c = str(r['seqnames'])
        if not c.startswith('chr'): c = 'chr' + c
        target_peaks.append((c, int(r['start']), int(r['end'])))
    print(f"      -> {len(sig)} aktive Refed-Peaks extrahiert.")

print(f" -> Insgesamt {len(target_peaks)} Peaks gesammelt.")

promoters = []
with open(base_bed_path) as f:
    for line in f:
        p = line.strip().split("\t")
        if len(p) >= 6:
            tss = int(p[1]) if p[5] == '+' else int(p[2])
            promoters.append((p[0], max(0, tss - 100), tss + 100, line))

matched = {
    g_line for g_chr, g_s, g_e, g_line in promoters
    if any(g_chr == p_chr and max(g_s, p_s) < min(g_e, p_e) for p_chr, p_s, p_e in target_peaks)
}

if not matched:
    print(" -> FEHLER: 0 Treffer nach Gen-Abgleich.")
    sys.exit(1)

print(f" -> EXAKT {len(matched)} Gene entsprechen nun dem finalen Kern-Set.")
with open(output_bed_path, "w") as out:
    out.writelines(sorted(matched))
EOF

# ==============================================================================
# 4. PLOTS ERZEUGEN
# ==============================================================================
shopt -s nullglob
for target in Creba Pol2; do
    rep_files=("${BIGWIG_DIR}/${target}_ChIP_"*"_4h_rep"*.bw)
    [ ${#rep_files[@]} -gt 0 ] || continue

    matrix_file="${OUTPUT_DIR}/matrix/${target}_states_4h_stringent_matrix.gz"
    output_pdf="${OUTPUT_DIR}/plots/combined_states_stringent/${target}_states_overlay_4h_sd.pdf"

    echo "Berechne Matrix für ${target} (4h Zustände, alle Replikate)..."
    computeMatrix scale-regions \
      -S "${rep_files[@]}" \
      -R "${FILTERED_BED}" \
      --beforeRegionStartLength 2000 \
      --regionBodyLength 3000 \
      --afterRegionStartLength 2000 \
      --skipZeros \
      -out "${matrix_file}" \
      -p 8

    echo "Generiere SD-Overlay-Plot für ${target} (4h)..."
    "$PYTHON_BIN" - "${matrix_file}" "${output_pdf}" "${target} Target States (4h) ±1 SD" << 'EOF'
import os, re, sys, gzip, json
from collections import defaultdict
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
num_samples = len(sample_labels) if sample_labels else 1
bins_per_sample = (df.shape[1] - meta_cols) // num_samples

groups = defaultdict(list)
for i in range(num_samples):
    name = os.path.basename(sample_labels[i]).replace(".bw", "")
    cond = re.sub(r'_rep[0-9]+$', '', re.sub(r'^(CrebA|Pol2)_ChIP_', '', name))
    s = meta_cols + i * bins_per_sample
    mat = df.iloc[:, s:s + bins_per_sample].apply(pd.to_numeric, errors='coerce').to_numpy(dtype=float)
    groups[cond].append(mat)

x = np.arange(bins_per_sample)
colors = ['#e41a1c', '#377eb8', '#4daf4a', '#ff7f00', '#984ea3']

plt.figure(figsize=(9, 6))
for idx, (cond, mats) in enumerate(groups.items()):
    stack = np.stack(mats, axis=0)  # (n_reps, n_genes, n_bins)
    n_reps = stack.shape[0]
    color = colors[idx % len(colors)]

    mean_profile = np.nanmean(np.nanmean(stack, axis=0), axis=0)
    sd_profile = np.nanmean(np.nanstd(stack, axis=0, ddof=1), axis=0) if n_reps > 1 else np.zeros_like(mean_profile)

    plt.plot(x, mean_profile, label=f"{cond} (n={n_reps})", color=color, linewidth=2)
    if n_reps > 1:
        plt.fill_between(x, mean_profile - sd_profile, mean_profile + sd_profile, color=color, alpha=0.22, linewidth=0)

tss_idx = int(bins_per_sample * (2000 / 7000))
tes_idx = int(bins_per_sample * (5000 / 7000))

plt.axvline(x=tss_idx, color='gray', linestyle='--', alpha=0.7, label='TSS')
plt.axvline(x=tes_idx, color='gray', linestyle=':', alpha=0.7, label='TES')
plt.xticks([0, tss_idx, tes_idx, bins_per_sample], ["-2.0 kb", "TSS", "TES", "+2.0 kb"])

plt.title(title, fontsize=14, fontweight='bold')
plt.xlabel("Genomic Region", fontsize=12)
plt.ylabel("Mean RPM Coverage", fontsize=12)
plt.legend(frameon=True, facecolor='white', edgecolor='none', fontsize=10)
plt.tight_layout()
plt.savefig(output_path, dpi=300)
plt.close()
print(f"SD-Overlay-Plot gespeichert unter: {output_path}")
EOF
done
shopt -u nullglob

echo "========================================================"
echo " Fertig! Nur SD-Plots erzeugt, alte Plots entfernt."
echo "========================================================"