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

# Filtere nur ChIP-Proben für das Peak-Calling heraus
if "Sample_Type" in metadata.columns:
    metadata["Sample_Type"] = metadata["Sample_Type"].astype(str).str.strip()
    CHIP_SAMPLES = metadata[metadata["Sample_Type"].str.upper() == "CHIP"]["Sample_ID"].tolist()
else:
    CHIP_SAMPLES = metadata["Sample_ID"].tolist()

metadata.set_index("Sample_ID", drop=False, inplace=True)
SAMPLES = metadata["Sample_ID"].tolist()

# 3. HILFSFUNKTIONEN (MÜSSEN IMMER VOR DEN REGELN STEHEN!)

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

def get_control_bam(wildcards):
    sample_id = wildcards.sample
    if "Input_Sample" in metadata.columns:
        control_id = metadata.loc[sample_id, "Input_Sample"]
        if pd.notna(control_id) and str(control_id).strip() != "" and str(control_id).upper() != "NAN":
            return f"results/bam/{control_id}.sorted.bam"
    raise ValueError(f"Keine gültige Input_Sample für Probe {sample_id} in Metadata.txt gefunden.")

# Liste aller Peak-Dateien für die jeweilige Gruppe
def get_group_peaks(wildcards):
    group_samples = chip_metadata[chip_metadata["Consensus_Group"] == wildcards.group]["Sample_ID"].tolist()
    return [f"results/peaks/narrowPeak/{sample}_peaks.narrowPeak" for sample in group_samples]

# 4. REGELN & DYNAMISCHE GRUPPEN
chip_metadata = metadata[metadata["Sample_Type"].str.upper() == "CHIP"].copy()

# Baut den Gruppennamen automatisch aus Target, Condition und (falls vorhanden) formatiertem Timepoint zusammen
time_suffix = "_" + chip_metadata["Timepoint"].astype(str).str.replace(r'h$', '', case=False, regex=True) + "h" if "Timepoint" in chip_metadata.columns else ""
chip_metadata["Consensus_Group"] = chip_metadata["Target"].astype(str) + "_" + chip_metadata["Condition"].astype(str) + time_suffix

CONSENSUS_GROUPS = chip_metadata["Consensus_Group"].unique().tolist()

rule all:
    input:
        expand("results/bigwig/{sample}.bw", sample=SAMPLES),
        expand("results/rds/{sample}_granges.rds", sample=SAMPLES),
        expand("results/peaks/narrowPeak/{sample}_peaks.narrowPeak", sample=CHIP_SAMPLES),
        expand("results/peaks/consensus/rds/{group}_consensus.rds", group=CONSENSUS_GROUPS),
        expand("results/peaks/consensus/bed/{group}_consensus.bed", group=CONSENSUS_GROUPS),
        expand("results/annotation/{group}_annotated_peaks.csv", group=CONSENSUS_GROUPS),
        expand("results/annotation/plots/{group}_anno_bar.pdf", group=CONSENSUS_GROUPS),
        expand("results/annotation/plots/{group}_dist_to_tss.pdf", group=CONSENSUS_GROUPS),
        expand("results/annotation/plots/{group}_dist_to_tss.pdf", group=CONSENSUS_GROUPS),
        # Universelle DiffBind Ausgaben
        "results/diffbind/differential_peaks.csv",
        "results/diffbind/plots/pca_plot.pdf",
        "results/diffbind/plots/ma_plot.pdf"

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
        frag_len = lambda wildcards: int(metadata.loc[wildcards.sample, 'Fragment_Le']) if "Fragment_Le" in metadata.columns and pd.notna(metadata.loc[wildcards.sample, 'Fragment_Le']) else 200
    log:
        "logs/r_coverage/{sample}.log"
    script:
        "scripts/calculate_coverage.R"

# 7. MACS2 Peak Calling für Drosophila (dm6) mit sortierten Unterordnern
rule macs2_callpeak:
    input:
        treatment = "results/bam/{sample}.sorted.bam",
        control = get_control_bam
    output:
        narrow = "results/peaks/narrowPeak/{sample}_peaks.narrowPeak",
        xls    = "results/peaks/xls/{sample}_peaks.xls",
        summits= "results/peaks/summits/{sample}_summits.bed"
    params:
        name = "{sample}",
        tmp_outdir = "results/peaks/tmp_{sample}",
        gsize = "dm",
        format_arg = lambda wildcards: "BAMPE" if metadata.loc[wildcards.sample, "Seq_Mode"] == "PE" else "BAM",
        extsize = lambda wildcards: f"--extsize {int(metadata.loc[wildcards.sample, 'Fragment_Le'])}" if metadata.loc[wildcards.sample, "Seq_Mode"] == "SE" and "Fragment_Le" in metadata.columns and pd.notna(metadata.loc[wildcards.sample, 'Fragment_Le']) else ""
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
            {params.extsize} \
            -q 0.01 \
            2> {log}
            
        mv {params.tmp_outdir}/{params.name}_peaks.narrowPeak {output.narrow}
        mv {params.tmp_outdir}/{params.name}_peaks.xls {output.xls}
        mv {params.tmp_outdir}/{params.name}_summits.bed {output.summits}
        
        rm -rf {params.tmp_outdir}
        """

# 8. KonsensPeaks 
rule make_consensus:
    input:
        get_group_peaks
    output:
        rds = "results/peaks/consensus/rds/{group}_consensus.rds",
        bed = "results/peaks/consensus/bed/{group}_consensus.bed"
    log:
        "logs/consensus/{group}.log"
    script:
        "scripts/make_consensus.R"

# 9. ChIPseeker Genom-Annotation und Visualisierung
rule chipseeker_annotation:
    input:
        rds = "results/peaks/consensus/rds/{group}_consensus.rds"
    output:
        csv = "results/annotation/{group}_annotated_peaks.csv",
        bar_pdf = "results/annotation/plots/{group}_anno_bar.pdf",
        tss_pdf = "results/annotation/plots/{group}_dist_to_tss.pdf"
    log:
        "logs/annotation/{group}.log"
    script:
        "scripts/run_chipseeker.R"

# 10. Differential Binding Analyse mit DiffBind (universell & korrigiert)
rule run_diffbind:
    input:
        bams = expand("results/bam/{sample}.sorted.bam", sample=CHIP_SAMPLES),
        peaks = expand("results/peaks/narrowPeak/{sample}_peaks.narrowPeak", sample=CHIP_SAMPLES)
    output:
        csv = "results/diffbind/differential_peaks.csv",
        pca = "results/diffbind/plots/pca_plot.pdf",
        ma  = "results/diffbind/plots/ma_plot.pdf"
    log:
        "logs/diffbind/diffbind.log"
    run:
        sample_ids_str = ",".join([f'"{s}"' for s in CHIP_SAMPLES])
        conditions_str = ",".join([f'"{metadata.loc[s, "Condition"]}"' for s in CHIP_SAMPLES])
        bams_str       = ",".join([f'"{x}"' for x in input.bams])
        peaks_str      = ",".join([f'"{x}"' for x in input.peaks])
        
        r_code = f"""
        library(DiffBind)
        
        samples <- data.frame(
          SampleID   = c({sample_ids_str}),
          Condition  = c({conditions_str}),
          bamReads   = c({bams_str}),
          Peaks      = c({peaks_str}),
          PeakCaller = "narrow"
        )

        db <- dba(sampleSheet = samples)
        db <- dba.count(db)
        
        # minMembers = 2 Replikate
        db <- dba.contrast(db, minMembers = 2)
        
        # dba.analyze braucht jetzt kein minMembers mehr
        db <- dba.analyze(db)

        res_table <- dba.report(db, th = 1)
        write.csv(as.data.frame(res_table), "{output.csv}", row.names = FALSE)
        
        pdf("{output.pca}"); dba.plotPCA(db, label = DBA_ID); dev.off()
        pdf("{output.ma}");  dba.plotMA(db); dev.off()
        """
        
        shell(f"Rscript -e '{r_code}' > {log} 2>&1")
