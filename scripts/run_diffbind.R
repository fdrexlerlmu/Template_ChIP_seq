options(error = function() {
  message("\n[FATAL ERROR STACK TRACE]")
  traceback(2)
  quit(status = 1)
})

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(GenomicRanges)
  library(DESeq2)
  library(Rsamtools)
  library(GenomeInfoDb)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("[ERROR] Zu wenig Argumente übergeben! Erwartet: 4, Erhalten: ", length(args))
}

meta_path       <- args[1]
output_pca      <- args[2]
reference_group <- args[3] 
target_name     <- args[4]

if (is.na(reference_group)) {
  reference_group <- "Fed_4h"
}

message("Verwende Target: ", target_name)
message("Verwende Referenzgruppe: ", reference_group)

# 1. Metadaten einlesen und bereinigen
meta <- read.delim(meta_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", check.names = FALSE)
if (ncol(meta) == 1) {
  meta <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
}
colnames(meta) <- gsub("[\r\n\ufeff]", "", trimws(colnames(meta)))
meta[] <- lapply(meta, function(x) if (is.character(x)) gsub("[\r\n\ufeff]", "", trimws(x)) else x)

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
col_tissue      <- find_col(meta, c("Tissue"))

meta_chip <- meta[tolower(meta[[col_sample_type]]) == "chip" & meta[[col_target]] == target_name, ]
meta_chip <- meta_chip[!duplicated(meta_chip[[col_sample_id]]), ]

if (nrow(meta_chip) == 0) {
  stop(paste("Keine Proben für Target", target_name, "in Metadata.txt gefunden!"))
}

clean_ids <- gsub("[^a-zA-Z0-9_-]", "", meta_chip[[col_sample_id]])
meta_chip$Group <- paste0(meta_chip[[col_condition]], "_", meta_chip[[col_time]])

bam_reads   <- file.path("results", "bam", paste0(clean_ids, ".sorted.bam"))
peaks_files <- file.path("results", "peaks", "narrowPeak", paste0(clean_ids, "_peaks.narrowPeak"))

for (f in c(bam_reads, peaks_files)) {
  if (!file.exists(f)) {
    stop(paste("[ERROR] Datei existiert nicht:", f))
  }
}

canonical_chroms <- c("2L", "2R", "3L", "3R", "4", "X", "Y")

# 2. Konsensus-Peakset via Base-R Tabellenimport erstellen
peak_list <- lapply(peaks_files, function(pf) {
  if (file.info(pf)$size == 0) {
    message("[WARNING] Peak-Datei ist leer: ", pf)
    return(GRanges())
  }
  df <- read.delim(pf, header = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(GRanges())
  
  chroms <- sub("^chr", "", as.character(df[, 1]))
  starts <- as.integer(df[, 2])
  ends <- as.integer(df[, 3])
  
  gr <- GRanges(seqnames = chroms, ranges = IRanges(start = starts + 1, end = ends))
  gr <- keepSeqlevels(gr, canonical_chroms[canonical_chroms %in% seqlevels(gr)], pruning.mode = "coarse")
  gr
})

unlist_peaks <- unlist(as(peak_list, "GRangesList"))
if (length(unlist_peaks) == 0) {
  stop("[ERROR] Alle Peak-Dateien sind leer oder enthalten keine Intervalle auf kanonischen Chromosomen für Target: ", target_name)
}

consensus_peaks <- reduce(unlist_peaks)

bfl <- BamFileList(bam_reads, yieldSize = 5000000)
seqinfo(consensus_peaks) <- seqinfo(bfl[[1]])
consensus_peaks <- keepSeqlevels(consensus_peaks, canonical_chroms[canonical_chroms %in% seqlevels(consensus_peaks)], pruning.mode = "coarse")

if (length(consensus_peaks) == 0) {
  stop("[ERROR] Konsensus-Peakset hat nach Seqinfo-Filterung 0 Intervalle.")
}

message("Konsensus-Peakset erstellt. Anzahl Intervalle: ", length(consensus_peaks))

# 3. Zählung im seriellen Modus direkt über BPPARAM (umgeht register() Bug)
message("Starte summarizeOverlaps (seriell)...")

se <- summarizeOverlaps(
  consensus_peaks, 
  bfl, 
  mode = "Union", 
  singleEnd = TRUE, 
  ignore.strand = TRUE,
  BPPARAM = SerialParam()
)

count_matrix <- assay(se)
colnames(count_matrix) <- clean_ids
count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

# 4. DESeq2 Analyse
col_data <- data.frame(
  Condition = factor(meta_chip$Group),
  row.names = clean_ids
)

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ Condition)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

target_out_dir  <- file.path("results/diffbind", target_name)
target_plot_dir <- file.path(target_out_dir, "plots")
dir.create(target_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(target_plot_dir, recursive = TRUE, showWarnings = FALSE)

# 5. PCA Plot
pdf(output_pca, width = 8, height = 6)
rld <- rlog(dds, blind = FALSE)
print(plotPCA(rld, intgroup = "Condition"))
dev.off()

# 6. Kontraste berechnen & exportieren
all_groups <- levels(col_data$Condition)
test_groups <- all_groups[all_groups != reference_group]

for (g in test_groups) {
  contrast_name <- paste0(g, "_vs_", reference_group)
  message("Berechne Kontrast: ", contrast_name)
  
  res_g <- results(dds, contrast = c("Condition", g, reference_group))
  res_df <- cbind(as.data.frame(consensus_peaks), as.data.frame(res_g))
  
  out_csv_path <- file.path(target_out_dir, paste0("differential_peaks_", contrast_name, ".csv"))
  write.csv(res_df, file = out_csv_path, row.names = FALSE)
  
  ma_path <- file.path(target_plot_dir, paste0("ma_plot_", contrast_name, ".pdf"))
  pdf(ma_path, width = 7, height = 6)
  plotMA(res_g, ylim = c(-4, 4), main = paste("MA Plot —", target_name, "|", contrast_name))
  abline(h = 0, col = "blue", lwd = 2)
  dev.off()
}

message("[✓] Analyse erfolgreich abgeschlossen für Target: ", target_name)