#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------------------
# 1. Pfade, Referenzen & Metadaten
# ------------------------------------------------------------------------------
RAW_DIR="rawdata"
REF_DIR="ref"
BAM_DIR="results/bam"
BW_DIR="results/bigwig"
RDS_DIR="results/rds"
METADATA="Metadata.txt"

# Referenzgenom (Drosophila melanogaster release-110 / dm6)
GENOME_URL="https://ftp.ensembl.org/pub/release-110/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz"
GENOME_FA="${REF_DIR}/genome.fa.gz"
INDEX_PREFIX="${REF_DIR}/bt2_index"

# CPU-Kerne ermitteln (mit 2 Kernen Puffer für flüssiges Arbeiten)
TOTAL_CORES=$(sysctl -n hw.ncpu)
THREADS=$(( TOTAL_CORES - 2 ))
if [ "$THREADS" -lt 1 ]; then THREADS=1; fi

# ------------------------------------------------------------------------------
# 2. Ordnerstruktur anlegen & Metadata-Prüfung
# ------------------------------------------------------------------------------
mkdir -p "${REF_DIR}" "${BAM_DIR}" "${BW_DIR}" "${RDS_DIR}"

if [ ! -f "${METADATA}" ]; then
    echo "[ERROR] Metadaten-Datei '${METADATA}' im Hauptverzeichnis nicht gefunden!"
    exit 1
fi

# ------------------------------------------------------------------------------
# 3. Referenzgenom herunterladen & Indizieren (einmalig)
# ------------------------------------------------------------------------------
if [ ! -f "${INDEX_PREFIX}.1.bt2" ] && [ ! -f "${INDEX_PREFIX}.1.bt2g" ]; then
    echo "[INFO] Lade Referenzgenom herunter..."
    curl -L -o "${GENOME_FA}" "${GENOME_URL}"
    
    echo "[INFO] Baue Bowtie2 Index..."
    bowtie2-build "${GENOME_FA}" "${INDEX_PREFIX}"
fi

# ------------------------------------------------------------------------------
# 4. Einlesen & Prozessieren über Metadata.txt
# ------------------------------------------------------------------------------
tr -d '\r' < "${METADATA}" | tail -n +2 | while IFS=$'\t' read -r SAMPLE_ID_RAW TARGET SAMPLE_TYPE INPUT_ID GENOTYPE CONDITION REPLICATE TISSUE ANTIBODY FRAG_LEN SEQ_MODE || [ -n "${SAMPLE_ID_RAW:-}" ]; do
    
    [ -z "${SAMPLE_ID_RAW:-}" ] && continue

    # Bereinige Sample_ID von evtl. mit eingegebenen Endungen
    SAMPLE_ID=$(echo "${SAMPLE_ID_RAW}" | sed -E 's/\.(fastq|txt)(\.gz)?$//')
    MODE=$(echo "${SEQ_MODE:-SE}" | tr '[:lower:]' '[:upper:]') # Standard: SE

    BAM_OUT="${BAM_DIR}/${SAMPLE_ID}.sorted.bam"
    BW_OUT="${BW_DIR}/${SAMPLE_ID}.bw"
    RDS_OUT="${RDS_DIR}/${SAMPLE_ID}_granges.rds"

    echo "======================================================================"
    echo "[INFO] Prozessiere Probe: ${SAMPLE_ID}"
    echo "[INFO] Mode: ${MODE} | Target: ${TARGET:-NA} | Condition: ${CONDITION:-NA}"
    echo "======================================================================"

    # 4a. Suche nach den Rohdaten (Single-End vs. Paired-End)
    FASTQ_1=""
    FASTQ_2=""

    if [ "${MODE}" == "PE" ]; then
        # Paired-End: Suche nach R1/R2 bzw. _1/_2
        for EXT in txt.gz fastq.gz; do
            if [ -f "${RAW_DIR}/${SAMPLE_ID}_R1.${EXT}" ] && [ -f "${RAW_DIR}/${SAMPLE_ID}_R2.${EXT}" ]; then
                FASTQ_1="${RAW_DIR}/${SAMPLE_ID}_R1.${EXT}"
                FASTQ_2="${RAW_DIR}/${SAMPLE_ID}_R2.${EXT}"
                break
            elif [ -f "${RAW_DIR}/${SAMPLE_ID}_1.${EXT}" ] && [ -f "${RAW_DIR}/${SAMPLE_ID}_2.${EXT}" ]; then
                FASTQ_1="${RAW_DIR}/${SAMPLE_ID}_1.${EXT}"
                FASTQ_2="${RAW_DIR}/${SAMPLE_ID}_2.${EXT}"
                break
            fi
        done

        if [ -z "${FASTQ_1}" ] || [ -z "${FASTQ_2}" ]; then
            echo "[WARN] Paired-End Dateien für '${SAMPLE_ID}' (R1/R2) in ${RAW_DIR}/ nicht gefunden. Überspringe..."
            continue
        fi
        echo "[INFO] PE-Dateien gefunden: ${FASTQ_1} & ${FASTQ_2}"

    else
        # Single-End
        for EXT in txt.gz fastq.gz; do
            if [ -f "${RAW_DIR}/${SAMPLE_ID}.${EXT}" ]; then
                FASTQ_1="${RAW_DIR}/${SAMPLE_ID}.${EXT}"
                break
            fi
        done

        if [ -z "${FASTQ_1}" ]; then
            echo "[WARN] Single-End Datei für '${SAMPLE_ID}' in ${RAW_DIR}/ nicht gefunden. Überspringe..."
            continue
        fi
        echo "[INFO] SE-Datei gefunden: ${FASTQ_1}"
    fi

    # 4b. Bowtie2 Alignment & BAM-Sortierung (mit MAPQ >= 20 Filter)
    if [ -f "${BAM_OUT}" ]; then
        echo "[SKIP] ${BAM_OUT} existiert bereits."
    else
        echo "[INFO] Starte Alignment (Bowtie2)..."
        if [ "${MODE}" == "PE" ]; then
            bowtie2 -x "${INDEX_PREFIX}" -1 "${FASTQ_1}" -2 "${FASTQ_2}" -p "${THREADS}" | \
                samtools view -bS -q 20 - | \
                samtools sort -@ 2 -m 1G -o "${BAM_OUT}" -
        else
            bowtie2 -x "${INDEX_PREFIX}" -U "${FASTQ_1}" -p "${THREADS}" | \
                samtools view -bS -q 20 - | \
                samtools sort -@ 2 -m 1G -o "${BAM_OUT}" -
        fi

        echo "[INFO] Indiziere BAM-Datei..."
        samtools index "${BAM_OUT}"
    fi

    # 4c. R-Block: Coverage, Mito-Filter, 1-bp RPM & GRanges-Export
    if [ -f "${BW_OUT}" ] && [ -f "${RDS_OUT}" ]; then
        echo "[SKIP] ${BW_OUT} und ${RDS_OUT} existieren bereits."
    else
        echo "[INFO] Berechne 1-bp RPM Coverage & speichere GRanges in R (${MODE}-Modus)..."
        
        Rscript --vanilla - <<EOF
suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(GenomicRanges)
  library(rtracklayer)
})

bam_file <- "${BAM_OUT}"
out_bw   <- "${BW_OUT}"
out_rds  <- "${RDS_OUT}"
seq_mode <- "${MODE}"
frag_len <- ${FRAG_LEN:-200}

param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE))

if (seq_mode == "PE") {
  # Paired-End: echte biologische Fragmente einlesen
  pairs <- readGAlignmentPairs(bam_file, param = param)
  gr_extended <- granges(pairs)
} else {
  # Single-End: Reads auf 5'-Start verankern und manuell verlängern
  reads <- readGAlignments(bam_file, param = param)
  gr_reads <- granges(reads)
  gr_extended <- resize(gr_reads, width = frag_len, fix = "start")
}

# Mitochondriale Reads aus dem Nenner filtern
gr_extended <- gr_extended[seqnames(gr_extended) != "dmel_mitochondrion_genome"]

# Chromosomenränder trimmen
gr_extended <- trim(gr_extended)

# 1. GRanges-Objekt (.rds) für spätere R-Analysen abspeichern
saveRDS(gr_extended, file = out_rds)

# 2. Exakte 1-bp Genom-Coverage berechnen & auf RPM normalisieren
cov <- coverage(gr_extended)
total_reads <- length(gr_extended)
cov_rpm <- cov * (1e6 / total_reads)

export.bw(cov_rpm, out_bw)
EOF
        
        echo "[SUCCESS] GRanges gespeichert: ${RDS_OUT}"
        echo "[SUCCESS] BigWig (1 bp, RPM) erstellt: ${BW_OUT}"
    fi

done

echo "======================================================================"
echo "[COMPLETE] Alle Proben aus ${METADATA} erfolgreich verarbeitet!"
echo "======================================================================"