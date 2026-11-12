import pandas as pd
import os

# 1. Konfiguration laden
configfile: "config.yml"

# 2. Metadaten einlesen & säubern (OHNE Index-Setzung, um flexibel nach Run_ID & Sample_ID zu filtern)
metadata_path = config.get("metadata", "Metadata.txt")
metadata = pd.read_csv(metadata_path, sep="\t").dropna(how="all")
metadata.columns = metadata.columns.str.strip()

# Flexibler Check für Seq_Mode
mode_col = None
for col in metadata.columns:
    clean_col = col.lower().replace("_", "").replace("-", "")
    if clean_col == "seqmode":
        mode_col = col
        break

if mode_col is None:
    raise KeyError(f"[ERROR] Spalte 'Seq_Mode' nicht in {metadata_path} gefunden!")

metadata.rename(columns={mode_col: "Seq_Mode"}, inplace=True)

# Säubern aller Text-Spalten
for col in ["Run_ID", "Sample_ID", "Seq_Mode", "Sample_Type", "Input_Sample_ID"]:
    if col in metadata.columns:
        metadata[col] = metadata[col].astype(str).str.strip()

metadata["Seq_Mode"] = metadata["Seq_Mode"].str.upper()

# 1. Eindeutige Liste aller technischen Läufe (Run_IDs) für Rule 5a
RUN_IDS = metadata["Run_ID"].unique().tolist()

# 2. Eindeutige Liste aller biologischen Proben (Sample_IDs) für Rule 5b, Coverage, etc.
SAMPLES = metadata["Sample_ID"].unique().tolist()

# 3. Eindeutige Liste aller ChIP-Proben (biologisch) für Peak Calling & DiffBind
chip_metadata = metadata[metadata["Sample_Type"].str.upper() == "CHIP"].copy()
CHIP_SAMPLES = chip_metadata["Sample_ID"].unique().tolist()


# 3. HILFSFUNKTIONEN

def get_fastq_inputs(wildcards):
    """Sucht die passenden FASTQ/TXT-Dateien basierend auf der Run_ID in rawdata/."""
    run_id = wildcards.run_id
    # Mode anhand der Run_ID aus dem DataFrame ermitteln
    mode_series = metadata.loc[metadata["Run_ID"] == run_id, "Seq_Mode"]
    if mode_series.empty:
        raise KeyError(f"Run_ID '{run_id}' nicht in Metadata.txt gefunden!")
    mode = mode_series.values[0]

    raw_dir = config.get("raw_dir", "rawdata")
    exts = ["fastq.gz", "fq.gz", "txt.gz"]
    
    if mode == "PE":
        patterns = [
            ("", ""),
            ("_R1", "_R2"),
            ("_R1_001", "_R2_001"),
            ("_1", "_2")
        ]
        for p1, p2 in patterns:
            for ext in exts:
                r1 = f"{raw_dir}/{run_id}{p1}.{ext}"
                r2 = f"{raw_dir}/{run_id}{p2}.{ext}"
                if os.path.exists(r1) and os.path.exists(r2):
                    return {"r1": r1, "r2": r2}
        raise FileNotFoundError(f"Keine Paired-End Dateien für Run_ID '{run_id}' in '{raw_dir}/' gefunden.")
    else:
        patterns = ["", "_R1", "_R1_001", "_1"]
        for p in patterns:
            for ext in exts:
                f = f"{raw_dir}/{run_id}{p}.{ext}"
                if os.path.exists(f):
                    return {"r1": f}
        raise FileNotFoundError(f"Keine Single-End Datei für Run_ID '{run_id}' in '{raw_dir}/' gefunden.")

def get_control_bam(wildcards):
    """Findet das zugehörige gemergte Input-BAM laut Metadata.txt."""
    sample_id = wildcards.sample
    input_series = metadata.loc[metadata["Sample_ID"] == sample_id, "Input_Sample_ID"]
    if input_series.empty:
        raise KeyError(f"Sample_ID '{sample_id}' nicht in Metadata.txt gefunden!")
    input_id = input_series.values[0]
    return f"results/bam/{input_id}.sorted.bam"


# 4. TARGET ALL RULE
rule all:
    input:
        # Final gemergte BAMs und BAIs (biologische Replikate)
        expand("results/bam/{sample}.sorted.bam", sample=SAMPLES),
        expand("results/bam/{sample}.sorted.bam.bai", sample=SAMPLES),
        
        # BigWig & GRanges Export
        expand("results/bigwig/{sample}.bw", sample=SAMPLES),
        expand("results/rds/{sample}_granges.rds", sample=SAMPLES),
        
        # MACS2 Peaks pro ChIP-Probe
        expand("results/peaks/narrowPeak/{sample}_peaks.narrowPeak", sample=CHIP_SAMPLES),
        
        # Differential Analysis Ausgaben
        "results/diffbind/differential_peaks.csv",
        "results/diffbind/plots/pca_plot.pdf",
        "results/diffbind/plots/ma_plot.pdf"


# 5a. Einzelne Technische Replikate alignen (landen strikt in results/bam/raw/)
rule align_tech_replicate:
    input:
        unpack(get_fastq_inputs)  # Nimmt wildcards.run_id
    output:
        bam = temp("results/bam/raw/{run_id}.sorted.bam")
    params:
        prefix = config.get("genome", {}).get("index_prefix", "genome/bowtie2/dm6") if isinstance(config.get("genome"), dict) else config.get("bowtie2_index", "genome/bowtie2/dm6"),
        mode   = lambda wildcards: metadata.loc[metadata["Run_ID"] == wildcards.run_id, "Seq_Mode"].values[0]
    threads: 6
    log:
        "logs/align/{run_id}.log"
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


# 5b. Technische Replikate pro biologischer Sample_ID zusammenführen (landen in results/bam/)
rule merge_tech_replicates:
    input:
        bams = lambda wildcards: [
            f"results/bam/raw/{run_id}.sorted.bam" 
            for run_id in metadata[metadata["Sample_ID"] == wildcards.sample]["Run_ID"]
        ]
    output:
        bam = "results/bam/{sample,\w+}.sorted.bam",
        bai = "results/bam/{sample,\w+}.sorted.bam.bai"
    threads: 4
    log:
        "logs/merge/{sample}.log"
    run:
        if not input.bams:
            raise ValueError(f"Keine Matching Run_IDs für Sample '{wildcards.sample}' in Metadata.txt gefunden!")

        # Falls nur 1 technisches Replikat existiert (_t1), direkt kopieren
        if len(input.bams) == 1:
            shell("cp {input.bams[0]} {output.bam}")
        else:
            bams_str = " ".join(input.bams)
            shell("samtools merge -@ {threads} -f {output.bam} " + bams_str + " 2> {log}")
        
        shell("samtools index {output.bam}")


# 6. R-Coverage calculation (BigWig)
rule r_coverage_export:
    input:
        bam = "results/bam/{sample}.sorted.bam",
        bai = "results/bam/{sample}.sorted.bam.bai"
    output:
        bw  = "results/bigwig/{sample}.bw",
        rds = "results/rds/{sample}_granges.rds"
    params:
        seq_mode = lambda wildcards: metadata.loc[metadata["Sample_ID"] == wildcards.sample, "Seq_Mode"].values[0],
        frag_len = 200
    log:
        "logs/r_coverage/{sample}.log"
    script:
        "scripts/calculate_coverage.R"


# 7. MACS2 Peak Calling
rule macs2_callpeak:
    input:
        treatment = "results/bam/{sample}.sorted.bam",
        control   = get_control_bam
    output:
        narrow  = "results/peaks/narrowPeak/{sample}_peaks.narrowPeak",
        xls     = "results/peaks/xls/{sample}_peaks.xls",
        summits = "results/peaks/summits/{sample}_summits.bed"
    params:
        name = "{sample}",
        tmp_outdir = "results/peaks/tmp_{sample}",
        gsize = "dm",
        format_arg = lambda wildcards: "BAMPE" if metadata.loc[metadata["Sample_ID"] == wildcards.sample, "Seq_Mode"].values[0] == "PE" else "BAM"
    log:
        "logs/macs2/{sample}.log"
    shell:
        """
        mkdir -p {params.tmp_outdir}
        
        macs2 callpeak \
            -t {input.treatment} \
            -c {input.control} \
            -f {params.format_arg} \
            -g {params.gsize} \
            -n {params.name} \
            --outdir {params.tmp_outdir} \
            -q 0.01 \
            2> {log}
            
        mv {params.tmp_outdir}/{params.name}_peaks.narrowPeak {output.narrow}
        mv {params.tmp_outdir}/{params.name}_peaks.xls {output.xls}
        mv {params.tmp_outdir}/{params.name}_summits.bed {output.summits}
        
        rm -rf {params.tmp_outdir}
        """

# 8. Differential Peak Analysis via Direct DESeq2 Pipeline
rule run_diffbind:
    input:
        meta  = metadata_path,
        bams  = expand("results/bam/{sample}.sorted.bam", sample=CHIP_SAMPLES),
        peaks = expand("results/peaks/narrowPeak/{sample}_peaks.narrowPeak", sample=CHIP_SAMPLES)
    output:
        csv = "results/diffbind/differential_peaks.csv",
        pca = "results/diffbind/plots/pca_plot.pdf",
        ma  = "results/diffbind/plots/ma_plot.pdf"
    log:
        "logs/diffbind/diffbind.log"
    shell:
        "Rscript --vanilla scripts/run_diffbind.R {input.meta} {output.csv} {output.pca} {output.ma} > {log} 2>&1"
