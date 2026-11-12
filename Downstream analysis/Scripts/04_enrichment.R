#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Dm.eg.db) # Annotation für Drosophila melanogaster
  library(cli)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
anno_csv <- args[1]
out_dir  <- args[2]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cli_h1("Downstream: Functional Enrichment Analysis (clusterProfiler)")

df <- read.csv(anno_csv)

# 1. Entrez-IDs & Symbols für signifikante Up-Zielgene extrahieren
up_genes <- df %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  pull(ENTREZID) %>%
  as.character() %>%
  na.omit() %>%
  unique()

cli_alert_info("Verarbeite {length(up_genes)} eindeutige Entrez-Gen-IDs für die Anreicherungsanalyse.")

if (length(up_genes) >= 5) {
  
  # 2. GO Biological Process (BP) Enrichment
  cli_progress_step("Führe GO Biological Process Enrichment aus...", spinner = TRUE)
  ego_bp <- enrichGO(
    gene          = up_genes,
    OrgDb         = org.Dm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  cli_progress_done()
  
  # 3. GO Cellular Component (CC) Enrichment (Zeigt u.a. ER, Golgi, Ribosom)
  cli_progress_step("Führe GO Cellular Component Enrichment aus...", spinner = TRUE)
  ego_cc <- enrichGO(
    gene          = up_genes,
    OrgDb         = org.Dm.eg.db,
    keyType       = "ENTREZID",
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  cli_progress_done()
  
  # 4. KEGG Pathway Enrichment (Drosophila melanogaster = "dme")
  cli_progress_step("Führe KEGG Pathway Enrichment aus...", spinner = TRUE)
  kegg_res <- enrichKEGG(
    gene         = up_genes,
    organism     = "dme",
    keyType      = "kegg",
    pvalueCutoff = 0.05
  )
  cli_progress_done()
  
  # 5. Ergebnisse als RDS speichern (ohne Komprimierung für schnelles I/O)
  saveRDS(ego_bp, file.path(out_dir, "enrich_go_bp.rds"), compress = FALSE)
  saveRDS(ego_cc, file.path(out_dir, "enrich_go_cc.rds"), compress = FALSE)
  saveRDS(kegg_res, file.path(out_dir, "enrich_kegg.rds"), compress = FALSE)
  
  cli_alert_success("Enrichment-Analyse erfolgreich abgeschlossen!")
  
} else {
  cli_alert_warning("Zu wenige Zielgene (unter 5) für eine statisch valide Anreicherungsanalyse.")
}

cli_h1("Fertig")