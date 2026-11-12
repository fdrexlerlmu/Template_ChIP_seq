if (file.exists("renv/activate.R")) source("renv/activate.R")

library(GenomicAlignments)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(GenomeInfoDb)

bam_file <- snakemake@input[["bam"]]
bw_file  <- snakemake@output[["bw"]]
rds_file <- snakemake@output[["rds"]]
seq_mode <- snakemake@params[["seq_mode"]]

# Fragmentlänge absichern: Falls NA/NULL -> Fallback auf 200
frag_len <- suppressWarnings(as.integer(snakemake@params[["frag_len"]]))
if (is.na(frag_len) || length(frag_len) == 0 || frag_len <= 0) {
  frag_len <- 200
}

bamfile <- BamFile(bam_file)

# 1. Alignment einlesen & Fragment-Verlängerung anwenden
if (seq_mode == "PE") {
  param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE))
  alignments <- readGAlignmentPairs(bamfile, param = param)
  gr <- granges(alignments, use.mcols = FALSE)
} else {
  param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE))
  alignments <- readGAlignments(bamfile, param = param)
  gr <- granges(alignments, use.mcols = FALSE)
  
  # Single-End Reads in 5' -> 3' Richtung auf frag_len (200 bp) verlängern
  gr <- resize(gr, width = frag_len, fix = "start")
}

# 2. Hauptchromosomen filtern & UCSC-Nomenklatur erzwingen (2L -> chr2L)
keep_seqs <- c("2L", "2R", "3L", "3R", "4", "X", "Y", "chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")
gr <- gr[seqnames(gr) %in% keep_seqs]
seqlevels(gr) <- keep_seqs[keep_seqs %in% seqlevels(gr)]
seqlevelsStyle(gr) <- "UCSC"

# 3. Coverage (RPM-normalisiert) berechnen
mapped_reads <- length(gr)
if (mapped_reads == 0) {
  stop(paste("Keine gültigen Reads in BAM-Datei gefunden:", bam_file))
}

cov <- coverage(gr)
rpm_cov <- cov * (1e6 / mapped_reads)

# 4. Exportieren
export.bw(rpm_cov, con = bw_file)
saveRDS(gr, file = rds_file)