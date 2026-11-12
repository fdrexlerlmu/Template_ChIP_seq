if (file.exists("renv/activate.R")) source("renv/activate.R")

library(DiffBind)
library(DESeq2)
library(Rsamtools)

# Argumente aus dem Snakemake-Aufruf einlesen
args <- commandArgs(trailingOnly = TRUE)
meta_path  <- args[1]
output_csv <- args[2]
output_pca <- args[3]
output_ma  <- args[4]

# 1. Metadaten einlesen & auf biologische ChIP-Samples reduzieren
meta <- read.delim(meta_path, stringsAsFactors = FALSE)
meta_chip <- meta[trimws(tolower(meta$Sample_Type)) == "chip", ]

# Eindeutige biologische Replikate herausfiltern (1 Zeile pro Sample_ID)
meta_chip <- meta_chip[!duplicated(trimws(meta_chip$Sample_ID)), ]
n_rows <- nrow(meta_chip)

samples <- data.frame(
  SampleID   = trimws(as.character(meta_chip$Sample_ID)),
  Tissue     = trimws(as.character(meta_chip$Tissue)),
  Factor     = trimws(as.character(meta_chip$Target)),
  Condition  = trimws(as.character(meta_chip$Condition)),
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
res <- results(dds)

dir.create("results/diffbind/plots", recursive = TRUE, showWarnings = FALSE)

# 8. CSV Ergebnisse schreiben
res_df <- cbind(as.data.frame(counts_gr), as.data.frame(res))
write.csv(res_df, file = output_csv, row.names = FALSE)

# 9. PCA Plot
pdf(output_pca, width = 8, height = 6)
rld <- rlog(dds, blind = FALSE)
plotPCA(rld, intgroup = "Condition")
dev.off()

# 10. MA Plot
pdf(output_ma, width = 7, height = 6)
plotMA(res, ylim = c(-4, 4), main = "Fed vs. Fasted (DESeq2 Manual LibSize)")
abline(h = 0, col = "blue", lwd = 2)
dev.off()