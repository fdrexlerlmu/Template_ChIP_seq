#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Dm.eg.db)
  library(enrichplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Fehler: Bitte übergebe den Plot-Ordner des Kontrasts und das Enrichment-Ausgabe-Verzeichnis!")
}

plots_dir <- args[1] # Entspricht z.B. results/plots/Fasted_4h_vs_Fed_4h/
out_dir   <- args[2] # Entspricht z.B. results/enrichment/Fasted_4h_vs_Fed_4h/

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

contrast_name <- basename(plots_dir)
message("--> Führe GO-Enrichment aus für Kontrast: ", contrast_name)

rds_path <- file.path(plots_dir, "gene_lists.rds")
if (!file.exists(rds_path)) {
  message("[-] Keine gene_lists.rds in ", plots_dir, " gefunden. Überspringe Enrichment.")
  quit(status = 0)
}

gene_lists <- readRDS(rds_path)

# Schleife für Upregulated & Downregulated
for (direction in names(gene_lists)) {
  genes <- gene_lists[[direction]]
  if (length(genes) < 5) {
    message("[-] Zu wenige Gene (< 5) für Richtung: ", direction)
    next
  }
  
  # Gene-Symbols in Entrez-IDs konvertieren
  gene_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dm.eg.db)
  if (nrow(gene_ids) == 0) {
    message("[-] Konnte keine Entrez-IDs mappen für Richtung: ", direction)
    next
  }
  
  entrez_vector <- unique(gene_ids$ENTREZID)
  
  # --- A. GO BP (Biological Process) ---
  ego_bp <- tryCatch({
    enrichGO(gene = entrez_vector, OrgDb = org.Dm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  }, error = function(e) NULL)
  
  if (!is.null(ego_bp) && nrow(ego_bp) > 0) {
    saveRDS(ego_bp, file.path(out_dir, paste0("enrich_go_bp_", direction, ".rds")))
    pdf(file.path(out_dir, paste0("dotplot_go_bp_", direction, ".pdf")), width = 9, height = 7)
    print(dotplot(ego_bp, showCategory = 20, title = paste("GO BP —", contrast_name, "(", direction, ")")))
    dev.off()
    message("[✓] GO BP Enrichment gespeichert für ", direction)
  }
  
  # --- B. GO CC (Cellular Component) ---
  ego_cc <- tryCatch({
    enrichGO(gene = entrez_vector, OrgDb = org.Dm.eg.db, keyType = "ENTREZID", ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  }, error = function(e) NULL)
  
  if (!is.null(ego_cc) && nrow(ego_cc) > 0) {
    saveRDS(ego_cc, file.path(out_dir, paste0("enrich_go_cc_", direction, ".rds")))
    pdf(file.path(out_dir, paste0("barplot_go_cc_", direction, ".pdf")), width = 9, height = 7)
    print(barplot(ego_cc, showCategory = 10, title = paste("GO CC —", contrast_name, "(", direction, ")")))
    dev.off()
    message("[✓] GO CC Enrichment gespeichert für ", direction)
  }
}

message("[✓] Enrichment-Analysen für Kontrast '", contrast_name, "' erfolgreich abgeschlossen.")