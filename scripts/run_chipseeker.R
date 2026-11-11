suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicRanges)
  library(rtracklayer)
  library(GenomeInfoDb)
  library(org.Dm.eg.db)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
})
  
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

# Inputs und Outputs aus dem S4 Snakemake-Objekt auslesen
rds_input   <- snakemake@input[["rds"]]
output_csv  <- snakemake@output[["csv"]]
output_bar  <- snakemake@output[["bar_pdf"]]
output_tss  <- snakemake@output[["tss_pdf"]]
group_name  <- snakemake@wildcards[["group"]]

# Ziel-Ordner sicherstellen
dir.create(dirname(output_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(output_bar), showWarnings = FALSE, recursive = TRUE)

cat(sprintf("===> Annoniere Konsens-Peaks für Gruppe: %s\n", group_name))

# 1. Konsens-GRanges einlesen
consensus_gr <- readRDS(rds_input)

# 2. Chromosomen-Namen an TxDb anpassen
seqlevelsStyle(consensus_gr) <- seqlevelsStyle(txdb)[1]

# 3. Peak-Annotation durchführen (TSS-Region -1000 bis +100 bp)
peakAnno <- annotatePeak(consensus_gr, 
                         tssRegion = c(-1000, 100),
                         TxDb = txdb, 
                         annoDb = "org.Dm.eg.db")

# Als Dataframe exportieren
anno_df <- as.data.frame(peakAnno)
write.csv(anno_df, file = output_csv, row.names = FALSE)
cat(sprintf("  [GESPEICHERT] CSV Tabelle: %s\n", output_csv))

# 4. Standard-Visualisierungen generieren und als PDF speichern
pdf(output_bar, width = 8, height = 6)
print(plotAnnoBar(peakAnno, title = paste0("Genomic Annotation: ", group_name)))
dev.off()

pdf(output_tss, width = 8, height = 6)
print(plotDistToTSS(peakAnno, title = paste0("Distribution relative to TSS: ", group_name)))
dev.off()

cat(sprintf("  [PLOTS FERTIG] Annotationen für %s erfolgreich erstellt.\n", group_name))