#!/usr/bin/env Rscript

# Automatisch ins Hauptverzeichnis wechseln, falls das Skript aus Downstream_analysis gestartet wurde
if (file.exists("../config.yml") && !file.exists("config.yml")) {
  setwd("..")
  message("[INFO] Wechsle automatisch ins Hauptverzeichnis: ", getwd())
}

suppressPackageStartupMessages({
  library(yaml)
  library(S4Vectors)
  library(IRanges)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(Rsamtools)
  library(DESeq2)
  library(GenomeInfoDb)
})

# 1. Konfiguration und Metadaten einlesen
config_path <- "config.yml"
if (!file.exists(config_path)) {
  stop("[ERROR] Konnte 'config.yml' nicht finden! Aktuelles Verzeichnis: ", getwd())
}

config <- yaml.load_file(config_path)

meta_path <- if (!is.null(config$metadata)) config$metadata else "Metadata.txt"
if (!file.exists(meta_path)) {
  stop(sprintf("[ERROR] Konnte Metadatendatei '%s' nicht finden!", meta_path))
}

message("[INFO] Lade Metadaten von: ", meta_path)
meta <- read.delim(meta_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", check.names = FALSE)
colnames(meta) <- gsub("[\r\n\ufeff]", "", trimws(colnames(meta)))
meta[] <- lapply(meta, function(x) if (is.character(x)) gsub("[\r\n\ufeff]", "", trimws(x)) else x)

# Hilfsfunktion zum Spaltenfinden
find_col <- function(df, candidates) {
  found <- candidates[candidates %in% colnames(df)]
  if (length(found) > 0) return(found[1])
  for (c in candidates) {
    match_idx <- grep(c, colnames(df), ignore.case = TRUE)
    if (length(match_idx) > 0) return(colnames(df)[match_idx[1]])
  }
  stop(sprintf("[ERROR] Spalte ähnlich zu '%s' nicht gefunden!", candidates[1]))
}

col_sample_id   <- find_col(meta, c("Sample_ID", "SampleID"))
col_target      <- find_col(meta, c("Target", "Factor"))
col_sample_type <- find_col(meta, c("Sample_Type", "SampleType"))
col_condition   <- find_col(meta, c("Condition"))
col_time        <- find_col(meta, c("Time"))

reference_group <- if (!is.null(config$reference_condition)) config$reference_condition else "Fasted_4h"

meta_chip <- meta[tolower(meta[[col_sample_type]]) == "chip", ]
meta_chip$Group <- paste0(meta_chip[[col_condition]], "_", meta_chip[[col_time]])

targets <- unique(meta_chip[[col_target]])
message("[INFO] Gefundene Targets für Analyse: ", paste(targets, collapse = ", "))

canonical_chroms <- c("2L", "2R", "3L", "3R", "4", "X", "Y")

# 2. Schleife über alle Targets
for (target_name in targets) {
  message("-----------------------------------------------------")
  message("[INFO] Starte Analyse für Target: ", target_name)
  
  sub_meta <- meta_chip[meta_chip[[col_target]] == target_name, ]
  sub_meta <- sub_meta[!duplicated(sub_meta[[col_sample_id]]), ]
  
  clean_ids <- gsub("[^a-zA-Z0-9_-]", "", sub_meta[[col_sample_id]])
  bam_reads   <- file.path("results", "bam", paste0(clean_ids, ".sorted.bam"))
  peaks_files <- file.path("results", "peaks", "narrowPeak", paste0(clean_ids, "_peaks.narrowPeak"))
  
  missing_bams <- bam_reads[!file.exists(bam_reads)]
  if (length(missing_bams) > 0) {
    warning("[WARNING] Überspringe Target ", target_name, ", fehlende BAMs: ", paste(missing_bams, collapse = ", "))
    next
  }
  
  peak_list <- lapply(peaks_files, function(pf) {
    if (!file.exists(pf) || file.info(pf)$size == 0) return(GRanges())
    df <- read.delim(pf, header = FALSE, stringsAsFactors = FALSE)
    if (nrow(df) == 0) return(GRanges())
    
    chroms <- sub("^chr", "", as.character(df[, 1]))
    starts <- as.integer(df[, 2])
    ends <- as.integer(df[, 3])
    
    gr <- GRanges(seqnames = chroms, ranges = IRanges(start = starts + 1, end = ends))
    keepSeqlevels(gr, canonical_chroms[canonical_chroms %in% seqlevels(gr)], pruning.mode = "coarse")
  })
  
  unlist_peaks <- unlist(as(peak_list, "GRangesList"))
  if (length(unlist_peaks) == 0) {
    warning("[WARNING] Keine Peaks für Target gefunden: ", target_name)
    next
  }
  
  consensus_peaks <- reduce(unlist_peaks)
  
  bfl <- BamFileList(bam_reads, yieldSize = 2000000)
  bam_si <- seqinfo(bfl[[1]])
  
  if (any(grepl("^chr", seqlevels(bam_si))) && !any(grepl("^chr", seqlevels(consensus_peaks)))) {
    seqlevels(consensus_peaks) <- paste0("chr", seqlevels(consensus_peaks))
    canonical_chroms_matched <- paste0("chr", canonical_chroms)
  } else {
    canonical_chroms_matched <- if (any(grepl("^chr", seqlevels(bam_si)))) paste0("chr", canonical_chroms) else canonical_chroms
  }
  
  seqinfo(consensus_peaks, pruning.mode = "coarse") <- bam_si[seqlevels(consensus_peaks)]
  consensus_peaks <- keepSeqlevels(consensus_peaks, canonical_chroms_matched[canonical_chroms_matched %in% seqlevels(consensus_peaks)], pruning.mode = "coarse")
  consensus_peaks <- trim(consensus_peaks)
  
  message("[INFO] Konsensus-Peakset Größe für ", target_name, ": ", length(consensus_peaks), " Intervalle")
  
  se <- summarizeOverlaps(
    consensus_peaks, 
    bfl, 
    mode = "Union", 
    singleEnd = TRUE, 
    ignore.strand = TRUE
  )
  
  count_matrix <- assay(se)
  colnames(count_matrix) <- clean_ids
  count_matrix <- round(count_matrix)
  storage.mode(count_matrix) <- "integer"
  
  condition_factor <- factor(sub_meta$Group)
  if (!(reference_group %in% levels(condition_factor))) {
    stop(sprintf("[ERROR] Referenzgruppe '%s' ist nicht in den Bedingungen vorhanden (%s).", 
                 reference_group, paste(levels(condition_factor), collapse = ", ")))
  }
  condition_factor <- relevel(condition_factor, ref = reference_group)
  
  col_data <- data.frame(
    Condition = condition_factor,
    row.names = clean_ids
  )
  
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ Condition)
  dds <- estimateSizeFactors(dds)
  
  dds <- tryCatch({
    DESeq(dds, fitType = "local")
  }, error = function(e) {
    message("[WARNING] Lokaler Fit fehlgeschlagen für ", target_name, ", versuche parametric Fit...")
    DESeq(dds, fitType = "parametric")
  })
  
  target_out_dir  <- file.path("results", "diffbind", target_name)
  target_plot_dir <- file.path(target_out_dir, "plots")
  dir.create(target_out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(target_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  pdf(file.path(target_plot_dir, paste0("pca_", target_name, ".pdf")), width = 8, height = 6)
  tryCatch({
    rld <- rlog(dds, blind = FALSE)
    print(plotPCA(rld, intgroup = "Condition") + ggplot2::ggtitle(paste("PCA — Target:", target_name)))
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("PCA Plot Fehler:\n", e$message), col = "red")
  })
  dev.off()
  
  all_groups <- levels(col_data$Condition)
  test_groups <- all_groups[all_groups != reference_group]
  
  for (g in test_groups) {
    contrast_name <- paste0(g, "_vs_", reference_group)
    message("[INFO] Berechne Kontrast: ", contrast_name)
    
    res_g <- results(dds, contrast = c("Condition", g, reference_group))
    res_df <- as.data.frame(res_g)
    res_df$Peak_ID <- paste0(seqnames(consensus_peaks), ":", start(consensus_peaks), "-", end(consensus_peaks))
    res_df <- cbind(Peak_ID = res_df$Peak_ID, as.data.frame(consensus_peaks), res_df[, !(colnames(res_df) %in% "Peak_ID")])
    
    out_csv <- file.path(target_out_dir, paste0("differential_peaks_", contrast_name, ".csv"))
    write.csv(res_df, file = out_csv, row.names = FALSE)
    
    pdf(file.path(target_plot_dir, paste0("ma_plot_", contrast_name, ".pdf")), width = 7, height = 6)
    plotMA(res_g, ylim = c(-4, 4), main = paste("MA Plot —", target_name, "|", contrast_name))
    abline(h = 0, col = "blue", lwd = 2)
    dev.off()
  }
}

message("[✓] Alle Analysen erfolgreich abgeschlossen!")