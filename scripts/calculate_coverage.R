log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(GenomicRanges)
  library(rtracklayer)
  library(Rsamtools)
})

bam_file <- snakemake@input[["bam"]]
out_bw   <- snakemake@output[["bw"]]
out_rds  <- snakemake@output[["rds"]]
seq_mode <- snakemake@params[["seq_mode"]]
frag_len <- as.integer(snakemake@params[["frag_len"]])

param <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE))

if (seq_mode == "PE") {
  # Paired-End: echte biologische Fragmente einlesen
  pairs <- readGAlignmentPairs(bam_file, param = param)
  gr_extended <- granges(pairs)
} else {
  # Single-End: Reads auf 5'-Start verankern und auf Fragmentlänge verlängern
  reads <- readGAlignments(bam_file, param = param)
  gr_reads <- granges(reads)
  gr_extended <- resize(gr_reads, width = frag_len, fix = "start")
}

# WICHTIG: Seqinfo (Chromosomenlängen) direkt aus dem BAM-Header übertragen,
# damit trim() und coverage() exakt wissen, wo die Chromosomenenden liegen.
seqinfo(gr_extended) <- seqinfo(BamFile(bam_file))

# Mitochondriale Reads (dmel_mitochondrion_genome) herausfiltern
gr_extended <- gr_extended[seqnames(gr_extended) != "dmel_mitochondrion_genome"]

# Chromosomenränder sauber trimmen (jetzt mit bekannten Chromosomenlängen)
gr_extended <- trim(gr_extended)

# 1. GRanges-Objekt (.rds) abspeichern
saveRDS(gr_extended, file = out_rds)

# 2. Exakte 1-bp Genom-Coverage berechnen & auf RPM normalisieren
cov <- coverage(gr_extended)
total_reads <- length(gr_extended)
cov_rpm <- cov * (1e6 / total_reads)

export.bw(cov_rpm, out_bw)