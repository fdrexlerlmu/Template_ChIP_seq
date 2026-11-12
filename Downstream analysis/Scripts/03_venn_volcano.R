#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
anno_csv <- args[1]
out_dir <- args[2]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(anno_csv)

# 1. Gen-Listen pro Gruppe extrahieren (padj & log2FoldChange)
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

# 2. Volcano Data Pre-processing
df <- df %>%
  mutate(
    Significance = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                              ~ "Not Sig"
    )
  )

saveRDS(df, file.path(out_dir, "volcano_data.rds"))

message("[✓] Venn- & Volcano-Datensätze vorbereitet.")