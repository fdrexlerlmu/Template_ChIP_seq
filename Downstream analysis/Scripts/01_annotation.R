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
diffbind_csv <- args[1]
out_dir <- args[2]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. DiffBind CSV einlesen
if (!file.exists(diffbind_csv)) {
  stop(paste("Datei nicht gefunden:", diffbind_csv))
}
df <- read.csv(diffbind_csv)

# 2. GRanges-Objekt erstellen
gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)

# 3. Chromosomen-Namen harmonisieren (UCSC Format erzwingen)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

# Falls die CSV-Chromosomen kein 'chr' haben, füge es hinzu / passe Styles an:
seqlevelsStyle(gr) <- "UCSC"
seqlevelsStyle(txdb) <- "UCSC"

# Kannonsche Chromosomen filtern (verhindert Fehler bei Unmapped/Scaffold-Contigs)
keep_seqs <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")
gr <- keepSeqlevels(gr, intersect(seqlevels(gr), keep_seqs), pruning.mode = "coarse")

# 4. Annotation mit ChIPseeker
peakAnno <- annotatePeak(
  gr, 
  TxDb = txdb, 
  annoDb = "org.Dm.eg.db",
  tssRegion = c(-3000, 3000),
  verbose = FALSE
)

# 5. Ergebnisse speichern
anno_df <- as.data.frame(peakAnno)
write.csv(anno_df, file.path(out_dir, "annotated_peaks.csv"), row.names = FALSE)
saveRDS(peakAnno, file.path(out_dir, "peakAnno.rds"))

message("[✓] Annotation erfolgreich abgeschlossen.")