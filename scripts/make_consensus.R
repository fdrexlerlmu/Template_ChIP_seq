library(GenomicRanges)
library(rtracklayer)

# 1. Getrennte Unterordner für Konsens-Formate anlegen
dir.create("results/peaks/granges", showWarnings = FALSE, recursive = TRUE)
dir.create("results/peaks/consensus/rds", showWarnings = FALSE, recursive = TRUE)
dir.create("results/peaks/consensus/bed", showWarnings = FALSE, recursive = TRUE)

# 2. Hilfsfunktion zum Einlesen von narrowPeak (0-based zu 1-based)
load_narrowPeak <- function(filepath) {
  df <- read.table(filepath, header = FALSE, sep = "\t")
  gr <- GRanges(
    seqnames = df$V1,
    ranges = IRanges(start = df$V2 + 1, end = df$V3),
    score = df$V7,
    signalValue = df$V7,
    pVal = df$V8,
    qVal = df$V9,
    summit = df$V10 + 1
  )
  return(gr)
}

# 3. Alle von Snakemake übergebenen Peak-Dateien einlesen
input_files <- snakemake@input
gr_list <- list()

for (filepath in input_files) {
  file_basename <- basename(filepath)
  sample_name <- sub("_peaks\\.narrowPeak$", "", file_basename)
  
  if (file.exists(filepath)) {
    gr <- load_narrowPeak(filepath)
    gr_list[[sample_name]] <- gr
    
    # Einzelne Proben-GRanges
    rds_sample_path <- file.path("results/peaks/granges", paste0(sample_name, "_granges.rds"))
    saveRDS(gr, file = rds_sample_path)
  }
}

# 4. Konsens-Peak-Set berechnen
output_rds <- snakemake@output[["rds"]]
output_bed <- snakemake@output[["bed"]]

if (length(gr_list) >= 2) {
  consensus <- reduce(subsetByOverlaps(gr_list[[1]], gr_list[[2]]))
  
  if (length(gr_list) > 2) {
    for (i in 3:length(gr_list)) {
      consensus <- reduce(subsetByOverlaps(consensus, gr_list[[i]]))
    }
  }
  
  # Abspeichern in die jeweiligen Unterordner
  saveRDS(consensus, file = output_rds)
  rtracklayer::export(consensus, con = output_bed, format = "BED")
  
  cat(sprintf("[ERFOLG] Konsens erstellt: %d robuste Peaks.\n", length(consensus)))
  
} else if (length(gr_list) == 1) {
  consensus <- gr_list[[1]]
  saveRDS(consensus, file = output_rds)
  rtracklayer::export(consensus, con = output_bed, format = "BED")
  cat("[HINWEIS] Nur 1 Probe in Gruppe gefunden – verwende Einzel-Peaks.\n")
}