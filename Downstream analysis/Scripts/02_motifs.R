#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
  library(memes)
  library(BiocParallel)
  library(cli)
  library(parallel)
})

# 1. CPU-Kerne dynamisch ermitteln (2 Kerne Puffer für macOS/RStudio)
total_cores <- parallel::detectCores(logical = FALSE)
num_workers <- max(1, total_cores - 2)

# 2. Parallelisierung für R & OpenMP/STREME setzen
register(MulticoreParam(workers = num_workers))
Sys.setenv(OMP_NUM_THREADS = num_workers) # Multithreading für STREME CLI

conda_bin <- file.path(Sys.getenv("CONDA_PREFIX"), "bin")
Sys.setenv(MEME_BIN = conda_bin)
options(memes.bin_path = conda_bin)

args <- commandArgs(trailingOnly = TRUE)
anno_csv <- args[1]
out_dir <- args[2]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cli_h1("Downstream: Motif Enrichment (STREME)")
cli_alert_info("System-Check: {total_cores} Kerne erkannt. Nutze dynamisch {num_workers} Kerne (via OMP_NUM_THREADS).")

# 3. Daten einlesen & filtern
cli_progress_step("Lade & sortiere annotated_peaks.csv...", spinner = TRUE)
df <- read.csv(anno_csv)

df_up <- df %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  arrange(padj) %>% 
  head(500)
cli_progress_done()

# 4. Sequenzen extrahieren
if (nrow(df_up) > 0) {
  cli_alert_info("Verarbeite {nrow(df_up)} Top-Peaks für die Motif-Suche.")
  
  cli_progress_step("Erstelle GRanges & reduziere Überlappungen...", spinner = TRUE)
  sig_up <- makeGRangesFromDataFrame(df_up, keep.extra.columns = TRUE)
  sig_up_clean <- reduce(sig_up)
  cli_progress_done()
  
  cli_progress_step("Extrahiere Fasta-Sequenzen aus dm6...", spinner = TRUE)
  seqs_up <- get_sequence(sig_up_clean, BSgenome.Dmelanogaster.UCSC.dm6)
  cli_progress_done()
  
  # 5. STREME-Aufruf (ohne 'cores'-Argument, nutzt OMP_NUM_THREADS)
  cli_h2("Starte STREME Motif-Suche ({num_workers} Cores via OpenMP)...")
  
  streme_res <- runStreme(
    seqs_up, 
    control = "shuffle", 
    outdir = file.path(out_dir, "streme_up"),
    nmotifs = 5,
    verbosity = 2
  )
  
  cli_progress_step("Speichere Ergebnisse...", spinner = TRUE)
  saveRDS(streme_res, file.path(out_dir, "streme_up.rds"), compress = FALSE)
  cli_progress_done()
  
  cli_alert_success("STREME Motif-Analyse erfolgreich abgeschlossen!")
  
} else {
  cli_alert_warning("Keine signifikanten Up-Peaks (padj < 0.05 & log2FoldChange > 0) gefunden.")
}

cli_h1("Fertig")