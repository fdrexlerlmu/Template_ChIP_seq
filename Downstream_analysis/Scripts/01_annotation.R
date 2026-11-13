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

# Ausgabeverzeichnis erstellen (entspricht jetzt z.B. results/annotation/Fasted_4h_vs_Fed_4h/)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Die spezifische DiffBind-Kontrast-CSV einlesen
if (!file.exists(diffbind_csv)) {
  stop(paste("Datei nicht gefunden:", diffbind_csv))
}
df <- read.csv(diffbind_csv)

if (nrow(df) == 0) {
  message("[-] Warnung: Die Datei ", diffbind_csv, " ist leer. Überspringe Annotation.")
  quit(status = 0)
}

message("[*] Annotiere Peaks für Kontrast-Datei: ", basename(diffbind_csv))

# 2. In GenomicRanges konvertieren
# Hinweis: Je nach DiffBind-Version heißen die Spalten seqnames/Chr, start, end. 
# Wir stellen sicher, dass Standard-Spalten genutzt werden.
if (!"seqnames" %in% colnames(df) && "Chr" %in% colnames(df)) {
  df <- df %>% rename(seqnames = Chr, start = Start, end = End)
}

gr_sub <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)

# 3. Chromosomen-Namen harmonisieren (UCSC Standard) & filtern
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
seqlevelsStyle(gr_sub) <- "UCSC"
seqlevelsStyle(txdb) <- "UCSC"

keep_seqs <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")
gr_sub <- keepSeqlevels(gr_sub, intersect(seqlevels(gr_sub), keep_seqs), pruning.mode = "coarse")

if (length(gr_sub) == 0) {
  message("[-] Keine Peaks nach Chromosomen-Filterung übrig.")
  quit(status = 0)
}

# 4. ChIPseeker Annotation durchführen
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

# 5. Diagramme für diesen spezifischen Kontrast generieren
contrast_name <- basename(out_dir)

pdf(file.path(out_dir, "peak_annotation_bar.pdf"), width = 8, height = 5)
print(plotAnnoBar(peakAnno) + labs(title = paste("Genomische Annotation —", contrast_name)))
dev.off()

pdf(file.path(out_dir, "dist_to_tss.pdf"), width = 8, height = 5)
print(plotDistToTSS(peakAnno, title = paste("Distanz zum TSS —", contrast_name)))
dev.off()

message("[✓] Annotation und Plots erfolgreich erstellt für: ", contrast_name)