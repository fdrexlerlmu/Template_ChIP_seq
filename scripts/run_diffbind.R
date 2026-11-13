if (file.exists("renv/activate.R")) source("renv/activate.R")

library(DiffBind)
library(DESeq2)
library(Rsamtools)

# Argumente aus dem Snakemake-Aufruf einlesen
args <- commandArgs(trailingOnly = TRUE)
meta_path       <- args[1]
output_pca      <- args[2]
reference_group <- args[3] # Wird aus der config.yml via Snakemake übergeben (z.B. "Fed_4h")

# Fallback, falls das Skript mal direkt in der Konsole aufgerufen wird
if (is.na(reference_group)) {
  reference_group <- "Fed_4h"
}

message("Verwende Referenzgruppe für Kontraste: ", reference_group)

# 1. Metadaten einlesen & auf biologische ChIP-Samples reduzieren
meta <- read.delim(meta_path, stringsAsFactors = FALSE)
meta_chip <- meta[trimws(tolower(meta$Sample_Type)) == "chip", ]

# Eindeutige biologische Replikate herausfiltern (1 Zeile pro Sample_ID)
meta_chip <- meta_chip[!duplicated(trimws(meta_chip$Sample_ID)), ]
n_rows <- nrow(meta_chip)

# WICHTIG: Condition und Time kombinieren, damit unterschiedliche Zeitpunkte 
# und Ernährungszustände (z.B. Fasted_4h vs Fed_4h) sauber getrennt werden!
meta_chip$Group <- paste(trimws(as.character(meta_chip$Condition)), 
                         trimws(as.character(meta_chip$Time)), sep = "_")

samples <- data.frame(
  SampleID   = trimws(as.character(meta_chip$Sample_ID)),
  Tissue     = trimws(as.character(meta_chip$Tissue)),
  Factor     = trimws(as.character(meta_chip$Target)),
  Condition  = trimws(as.character(meta_chip$Group)), # Nutzt die kombinierte Gruppe
  Replicate  = seq_len(n_rows),
  bamReads   = paste0("results/bam/", trimws(meta_chip$Sample_ID), ".sorted.bam"),
  Peaks      = paste0("results/peaks/narrowPeak/", trimws(meta_chip$Sample_ID), "_peaks.narrowPeak"),
  PeakCaller = "narrow",
  bamControl = paste0("results/bam/", trimws(meta_chip$Input_Sample_ID), ".sorted.bam"),
  stringsAsFactors = FALSE
)

# 2. DiffBind Zählen mit expliziten Roh-Read-Counts (DBA_SCORE_READS)
db <- dba(sampleSheet = samples)
db <- dba.count(db, minOverlap = 2, summits = FALSE, score = DBA_SCORE_READS)

# 3. Count Matrix & Genomic Ranges extrahieren
counts_gr <- dba.peakset(db, bRetrieve = TRUE)
count_matrix <- as.matrix(mcols(counts_gr))
count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

colnames(count_matrix) <- samples$SampleID

# 4. Total Mapped Reads aus den ChIP-BAMs ermitteln
get_mapped_reads <- function(bam_path) {
  sum(idxstatsBam(bam_path)$mapped)
}
lib_sizes <- sapply(samples$bamReads, get_mapped_reads)

# 5. DESeq2 Objekt manuell konstruieren
col_data <- data.frame(
  Condition = factor(samples$Condition),
  row.names = samples$SampleID
)
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ Condition)

# 6. Bibliotheksgrößen als Size Factors zuweisen
sizeFactors(dds) <- lib_sizes / mean(lib_sizes)

# 7. DESeq2 Differential Binding rechnen
dds <- DESeq(dds)

# Verzeichnisse für die Ergebnisse erstellen
dir.create("results/diffbind", recursive = TRUE, showWarnings = FALSE)
dir.create("results/diffbind/plots", recursive = TRUE, showWarnings = FALSE)

# 8. Globaler PCA Plot (über alle Proben hinweg)
pdf(output_pca, width = 8, height = 6)
rld <- rlog(dds, blind = FALSE)
print(plotPCA(rld, intgroup = "Condition"))
dev.off()

# 9. Dynamische Kontraste gegen die Referenz berechnen (Schleife für alle Testgruppen)
all_groups <- levels(col_data$Condition)
test_groups <- all_groups[all_groups != reference_group]

for (g in test_groups) {
  contrast_name <- paste0(g, "_vs_", reference_group)
  message("Berechne Kontrast: ", contrast_name)
  
  # Spezifischen Kontrast über DESeq2 extrahieren
  res_g <- results(dds, contrast = c("Condition", g, reference_group))
  
  # Zusammenführen mit den Genomic Ranges
  res_df <- cbind(as.data.frame(counts_gr), as.data.frame(res_g))
  
  # Eigene CSV-Datei pro Kontrast speichern (wird von der Downstream-Pipeline eingelesen)
  out_csv_path <- file.path("results/diffbind", paste0("differential_peaks_", contrast_name, ".csv"))
  write.csv(res_df, file = out_csv_path, row.names = FALSE)
  
  # Individuellen MA-Plot für diesen spezifischen Kontrast speichern
  ma_path <- file.path("results/diffbind/plots", paste0("ma_plot_", contrast_name, ".pdf"))
  pdf(ma_path, width = 7, height = 6)
  plotMA(res_g, ylim = c(-4, 4), main = paste("MA Plot —", contrast_name))
  abline(h = 0, col = "blue", lwd = 2)
  dev.off()
}

message("[✓] DiffBind & DESeq2 Analysen für alle Kontraste erfolgreich abgeschlossen.")