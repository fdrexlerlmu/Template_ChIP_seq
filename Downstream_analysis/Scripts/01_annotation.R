#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(ChIPseeker)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
  library(org.Dm.eg.db)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Fehler: Bitte übergebe den Pfad zur CSV-Datei und das Ausgabeverzeichnis als Argumente!")
}

diffbind_csv <- args[1]
out_dir <- args[2]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(diffbind_csv)) {
  stop(paste("Datei nicht gefunden:", diffbind_csv))
}
df <- read.csv(diffbind_csv)

if (nrow(df) == 0) {
  message("[-] Warnung: Die Datei ", diffbind_csv, " ist leer. Überspringe Annotation.")
  quit(status = 0)
}

message("[*] Annotiere Peaks für Kontrast-Datei: ", basename(diffbind_csv))

# 1. Spaltennamen harmonisieren (DiffBind Standard)
if (!"seqnames" %in% colnames(df) && "Chr" %in% colnames(df)) {
  df <- df %>% rename(seqnames = Chr, start = Start, end = End)
}

# 2. WICHTIG: Nur signifikante, differentielle Peaks filtern!
if ("FDR" %in% colnames(df)) {
  df_sig <- df %>% filter(FDR < 0.05)
} else if ("padj" %in% colnames(df)) {
  df_sig <- df %>% filter(padj < 0.05)
} else if ("p.value" %in% colnames(df)) {
  df_sig <- df %>% filter(p.value < 0.05)
} else if ("pvalue" %in% colnames(df)) {
  df_sig <- df %>% filter(pvalue < 0.05)
} else {
  message("[!] Keine klassische Statistik-Spalte gefunden. Nutze alle Zeilen.")
  df_sig <- df
}

# 3. Fallback, falls der FDR < 0.05 Filter zu wenige Zeilen liefert (< 20)
if (nrow(df_sig) < 20) {
  message("[-] Warnung: Nur ", nrow(df_sig), " signifikante Peaks gefunden (< 20). Nutze Top 500 nach Fold Change als Fallback.")
  
  # Suche primär nach LogFC, um kontrast-spezifische Top-Regulierungen zu erhalten
  lfc_col <- if ("log2FoldChange" %in% colnames(df)) "log2FoldChange" else if ("Fold" %in% colnames(df)) "Fold" else NULL
  
  if (!is.null(lfc_col)) {
    # Sortiere nach absolutem Fold Change absteigend (stärkste Up- und Down-Regulierung)
    df_sig <- df %>% arrange(desc(abs(.data[[lfc_col]]))) %>% head(500)
  } else {
    # Absoluter Notfall-Fallback, falls gar kein Fold Change da ist
    df_sig <- df %>% head(500)
  }
}

message("[*] Anzahl signifikanter Peaks für Annotation: ", nrow(df_sig))

# 4. In GenomicRanges konvertieren
gr_sub <- makeGRangesFromDataFrame(df_sig, keep.extra.columns = TRUE)

# 5. Chromosomen-Namen harmonisieren (UCSC Standard) & filtern
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
seqlevelsStyle(gr_sub) <- "UCSC"
seqlevelsStyle(txdb) <- "UCSC"

keep_seqs <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")
gr_sub <- keepSeqlevels(gr_sub, intersect(seqlevels(gr_sub), keep_seqs), pruning.mode = "coarse")

if (length(gr_sub) == 0) {
  message("[-] Keine Peaks nach Chromosomen-Filterung übrig.")
  quit(status = 0)
}

# 6. ChIPseeker Annotation durchführen
peakAnno <- annotatePeak(
  gr_sub, 
  TxDb = txdb, 
  annoDb = "org.Dm.eg.db",
  tssRegion = c(-3000, 3000),
  verbose = FALSE
)

# Tabellen direkt im kontrast-spezifischen Ordner speichern
anno_df <- as.data.frame(peakAnno)
write.csv(anno_df, file.path(out_dir, "annotated_peaks.csv"), row.names = FALSE)
saveRDS(peakAnno, file.path(out_dir, "peakAnno.rds"))

# 7. Diagramme für diesen spezifischen Kontrast generieren
contrast_name <- basename(out_dir)

pdf(file.path(out_dir, "peak_annotation_bar.pdf"), width = 8, height = 5)
print(plotAnnoBar(peakAnno) + labs(title = paste("Genomische Annotation (Sig. Peaks) —", contrast_name)))
dev.off()

pdf(file.path(out_dir, "dist_to_tss.pdf"), width = 8, height = 5)
print(plotDistToTSS(peakAnno, title = paste("Distanz zum TSS —", contrast_name)))
dev.off()

message("[✓] Annotation und Plots erfolgreich erstellt für: ", contrast_name)