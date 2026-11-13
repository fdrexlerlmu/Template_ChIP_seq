#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Fehler: Bitte übergebe den Annotations-Ordner des Kontrasts und das Ausgabe-Verzeichnis!")
}

anno_dir <- args[1] # Entspricht z.B. results/annotation/Fasted_4h_vs_Fed_4h/
out_dir  <- args[2] # Entspricht z.B. results/plots/Fasted_4h_vs_Fed_4h/

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

contrast_name <- basename(anno_dir)
message("--> Bereite Venn & Volcano Daten vor für Kontrast: ", contrast_name)

anno_csv <- file.path(anno_dir, "annotated_peaks.csv")
if (!file.exists(anno_csv)) {
  message("[-] Keine annotated_peaks.csv in ", anno_dir, " gefunden. Überspringe.")
  quit(status = 0)
}

df <- read.csv(anno_csv)

# 1. Gen-Listen für Venn (Upregulated vs Downregulated im Kontrast)
genes_up <- df %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  pull(SYMBOL) %>% 
  na.omit() %>% 
  unique()

genes_down <- df %>% 
  filter(padj < 0.05 & log2FoldChange < 0) %>% 
  pull(SYMBOL) %>% 
  na.omit() %>% 
  unique()

gene_lists <- list(
  Upregulated = genes_up,
  Downregulated = genes_down
)

saveRDS(gene_lists, file.path(out_dir, "gene_lists.rds"))

# 2. Volcano-Daten aufbereiten
df_volcano <- df %>%
  mutate(
    Significance = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                              ~ "Not Sig"
    )
  )

saveRDS(df_volcano, file.path(out_dir, "volcano_data.rds"))
message("[✓] Venn- & Volcano-Daten für Kontrast '", contrast_name, "' erfolgreich gespeichert.")