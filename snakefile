import pandas as pd
import os

# 1. Konfiguration laden
configfile: "config.yml"

# 2. Metadaten einlesen & säubern
metadata_path = config.get("metadata", "Metadata.txt")
metadata = pd.read_csv(metadata_path, sep="\t").dropna(how="all")
metadata.columns = metadata.columns.str.strip()

# Flexibler Check für Seq_Mode
mode_col = None
for col in metadata.columns:
    if col.lower().replace("_", "") == "seqmode":
        mode_col = col
        break

if mode_col is None:
    raise KeyError(
        f"[ERROR] Spalte 'Seq_Mode' nicht in {metadata_path} gefunden! "
        f"Vorhandene Spalten sind: {list(metadata.columns)}"
    )

metadata.rename(columns={mode_col: "Seq_Mode"}, inplace=True)

metadata["Sample_ID"] = metadata["Sample_ID"].astype(str).str.strip()
metadata["Sample_ID"] = metadata["Sample_ID"].str.replace(r"\.(fastq|txt)(\.gz)?$", "", regex=True)
metadata["Seq_Mode"] = metadata["Seq_Mode"].astype(str).str.strip().str.upper()

metadata.set_index("Sample_ID", drop=False, inplace=True)
SAMPLES = metadata["Sample_ID"].tolist()

# 3. Hilfsfunktion ZUERST definieren (vor allen Regeln!)
def get_fastq_inputs(wildcards):
    sample_id = wildcards.sample
    mode = metadata.loc[sample_id, "Seq_Mode"]
    raw_dir = config.get("raw_dir", "rawdata")
    
    if mode == "PE":
        for ext in ["txt.gz", "fastq.gz"]:
            r1, r2 = f"{raw_dir}/{sample_id}_R1.{ext}", f"{raw_dir}/{sample_id}_R2.{ext}"
            if os.path.exists(r1) and os.path.exists(r2):
                return {"r1": r1, "r2": r2}
            r1_alt, r2_alt = f"{raw_dir}/{sample_id}_1.{ext}", f"{raw_dir}/{sample_id}_2.{ext}"
            if os.path.exists(r1_alt) and os.path.exists(r2_alt):
                return {"r1": r1_alt, "r2": r2_alt}
        raise FileNotFoundError(f"Keine Paired-End Dateien für {sample_id} in {raw_dir}/ gefunden.")
    else:
        for ext in ["txt.gz", "fastq.gz"]:
            f = f"{raw_dir}/{sample_id}.{ext}"
            if os.path.exists(f):
                return {"r1": f}
        raise FileNotFoundError(f"Keine Single-End Datei für {sample_id} in {raw_dir}/ gefunden.")

# 4. Zielregel (rule all)
rule all:
    input:
        expand("results/bigwig/{sample}.bw", sample=SAMPLES),
        expand("results/rds/{sample}_granges.rds", sample=SAMPLES)

# 5. Alignment & BAM-Sortierung (Bowtie2 + Samtools, MAPQ >= 20)
rule align_and_sort:
    input:
        unpack(get_fastq_inputs)
    output:
        bam = "results/bam/{sample}.sorted.bam",
        bai = "results/bam/{sample}.sorted.bam.bai"
    params:
        prefix = config["genome"]["index_prefix"],
        mode = lambda wildcards: metadata.loc[wildcards.sample, "Seq_Mode"]
    threads: 6
    log:
        "logs/align/{sample}.log"
    run:
        if params.mode == "PE":
            shell(
                "bowtie2 -x {params.prefix} -1 {input.r1} -2 {input.r2} -p {threads} 2> {log} | "
                "samtools view -bS -q 20 - | "
                "samtools sort -@ 2 -m 1G -o {output.bam} -"
            )
        else:
            shell(
                "bowtie2 -x {params.prefix} -U {input.r1} -p {threads} 2> {log} | "
                "samtools view -bS -q 20 - | "
                "samtools sort -@ 2 -m 1G -o {output.bam} -"
            )
        shell(
            "samtools index {output.bam}"
        )

# 6. R-Coverage calculation & GRanges Export
rule r_coverage_export:
    input:
        bam = "results/bam/{sample}.sorted.bam",
        bai = "results/bam/{sample}.sorted.bam.bai"
    output:
        bw = "results/bigwig/{sample}.bw",
        rds = "results/rds/{sample}_granges.rds"
    params:
        seq_mode = lambda wildcards: metadata.loc[wildcards.sample, "Seq_Mode"],
        frag_len = lambda wildcards: int(metadata.loc[wildcards.sample, "Fragment_Length"]) if "Fragment_Length" in metadata.columns and pd.notna(metadata.loc[wildcards.sample, "Fragment_Length"]) else 200
    log:
        "logs/r_coverage/{sample}.log"
    script:
        "scripts/calculate_coverage.R"
