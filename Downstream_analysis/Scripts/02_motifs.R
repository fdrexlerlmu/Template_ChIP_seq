#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
  library(memes)
  library(BiocParallel)
  library(cli)
  library(parallel)
  library(universalmotif)
  library(ggseqlogo)
})

# 1. CPU-Kerne dynamisch ermitteln (2 Kerne Puffer für macOS/RStudio)
total_cores <- parallel::detectCores(logical = FALSE)
num_workers <- max(1, total_cores - 2)

# 2. Parallelisierung für R & OpenMP/STREME setzen
register(MulticoreParam(workers = num_workers))
Sys.setenv(OMP_NUM_THREADS = num_workers)

conda_bin <- file.path(Sys.getenv("CONDA_PREFIX"), "bin")
Sys.setenv(MEME_BIN = conda_bin)
options(memes.bin_path = conda_bin)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Fehler: Bitte übergebe den Annotations-Ordner und den Motiv-Ausgabeordner als Argumente!")
}

anno_dir <- args[1]
out_dir  <- args[2]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

contrast_name <- basename(anno_dir)
cli_h1("Downstream: Motif Enrichment (STREME) für Kontrast: {contrast_name}")
cli_alert_info("System-Check: {total_cores} Kerne erkannt. Nutze dynamisch {num_workers} Kerne.")

specific_anno_csv <- file.path(anno_dir, "annotated_peaks.csv")

if (!file.exists(specific_anno_csv)) {
  cli_alert_warning("Keine annotierten Peaks unter {specific_anno_csv} gefunden. Überspringe STREME.")
  saveRDS(NULL, file.path(out_dir, "streme_up.rds"), compress = FALSE)
  quit(status = 0)
}

# 3. Daten einlesen & dynamische Spalten-Erkennung
cli_progress_step("Lade & sortiere annotierte Peaks...", spinner = TRUE)
df <- read.csv(specific_anno_csv)

stat_col <- if ("padj" %in% colnames(df)) "padj" else if ("FDR" %in% colnames(df)) "FDR" else "pvalue"
lfc_col  <- if ("log2FoldChange" %in% colnames(df)) "log2FoldChange" else "Fold"

if (!all(c(stat_col, lfc_col) %in% colnames(df))) {
  cli_alert_warning("Erforderliche Statistik- oder LogFC-Spalten in {specific_anno_csv} nicht gefunden.")
  saveRDS(NULL, file.path(out_dir, "streme_up.rds"), compress = FALSE)
  quit(status = 0)
}

df_up <- df %>%  
  filter(.data[[lfc_col]] > 0) %>%  
  arrange(.data[[stat_col]]) %>%  
  head(300)

cli_progress_done()

if (nrow(df_up) >= 5) {
  cli_alert_info("Verarbeite {nrow(df_up)} Top-Peaks für die Motif-Suche.")
  
  cli_progress_step("Erstelle GRanges & reduziere Überlappungen...", spinner = TRUE)
  sig_up <- makeGRangesFromDataFrame(df_up, keep.extra.columns = TRUE)
  sig_up_clean <- reduce(sig_up)
  cli_progress_done()
  
  cli_progress_step("Extrahiere Fasta-Sequenzen aus dm6...", spinner = TRUE)
  seqs_up <- get_sequence(sig_up_clean, BSgenome.Dmelanogaster.UCSC.dm6)
  cli_progress_done()
  
  streme_out_dir <- file.path(out_dir, "streme_up")
  dir.create(streme_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 4. STREME-Aufruf
  cli_alert_info("Starte STREME Motif-Suche...")
  
  streme_res <- tryCatch({
    runStreme(
      seqs_up, 
      control = "shuffle", 
      outdir = streme_out_dir,
      nmotifs = 5,
      verbosity = 2
    )
  }, error = function(e) {
    cli_alert_warning(paste("STREME-Lauf fehlgeschlagen:", e$message))
    return(NULL)
  })
  
  # NEU & ROBUST: Generiere Logo-PNGs direkt aus der von STREME geschriebenen streme.txt
  cli_progress_step("Generiere Sequenzlogo-PNGs...", spinner = TRUE)
  logo_dir <- file.path(streme_out_dir, "logos")
  dir.create(logo_dir, showWarnings = FALSE, recursive = TRUE)
  
  streme_txt <- file.path(streme_out_dir, "streme.txt")
  if (file.exists(streme_txt)) {
    motifs <- tryCatch(universalmotif::read_meme(streme_txt), error = function(e) NULL)
    
    if (!is.null(motifs)) {
      # Konvertiere in eine Liste (falls nur 1 Motiv gefunden wurde)
      motifs_list <- if (inherits(motifs, "universalmotif")) list(motifs) else motifs
      
      if (length(motifs_list) > 0) {
        for (idx in seq_along(motifs_list)) {
          m <- motifs_list[[idx]]
          if (!is.null(m@motif)) {
            p <- ggseqlogo(m@motif, method = 'bits') + 
              theme_minimal() + 
              labs(title = paste0("Motif ", idx))
            
            ggsave(
              filename = file.path(logo_dir, paste0("motif_", idx, ".png")), 
              plot = p, width = 5, height = 2.5, dpi = 300, bg = "white"
            )
          }
        }
      }
    }
  }
  cli_progress_done()
  
  cli_progress_step("Speichere Ergebnisse...", spinner = TRUE)
  saveRDS(streme_res, file.path(out_dir, "streme_up.rds"), compress = FALSE)
  cli_progress_done()
  
  cli_alert_success("STREME Motif-Analyse für {contrast_name} erfolgreich abgeschlossen!")
  
} else {
  cli_alert_warning("Zu wenige signifikante Up-Peaks ({nrow(df_up)} < 5) für {contrast_name}. STREME wird übersprungen.")
  saveRDS(NULL, file.path(out_dir, "streme_up.rds"), compress = FALSE)
}

cli_h1("Motif-Analyse für {contrast_name} abgeschlossen!")