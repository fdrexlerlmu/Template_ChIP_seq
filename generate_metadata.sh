#!/usr/bin/env bash
set -euo pipefail

RAW_DIR="rawdata"
OUTPUT_FILE="Metadata.txt"

if [ ! -d "$RAW_DIR" ]; then
    echo "[ERROR] Ordner '$RAW_DIR' nicht gefunden!" >&2
    exit 1
fi

tmp_list=$(mktemp)

for f in "$RAW_DIR"/*.gz; do
    [ -e "$f" ] || continue
    filename=$(basename "$f" | tr -d '\r')
    file_id=$(echo "$filename" | sed -E 's/(\.(fastq|fq|txt|gz))+$$//i')
    echo "$file_id" >> "$tmp_list"
done

# Header um die Spalte 'Time' ergänzt
echo -e "Run_ID\tSample_ID\tTarget\tSample_Type\tInput_Sample_ID\tGenotype\tCondition\tTime\tBio_Replicate\tTech_Replicate\tTissue\tSeq_Mode\tFragment_Length" > "$OUTPUT_FILE"

while read -r file_id; do
    [ -z "$file_id" ] && continue
    file_id=$(echo "$file_id" | tr -d '\r')

    num_fields=$(echo "$file_id" | awk -F'_' '{print NF}')

    if [ "$num_fields" -eq 6 ]; then
        target_raw=$(echo "$file_id" | awk -F'_' '{print $1}')
        stype_raw=$(echo "$file_id" | awk -F'_' '{print $2}')
        cond_raw=$(echo "$file_id" | awk -F'_' '{print $3}')
        time_raw=$(echo "$file_id" | awk -F'_' '{print $4}')
        bio_raw=$(echo "$file_id" | awk -F'_' '{print $5}')
        tech_raw=$(echo "$file_id" | awk -F'_' '{print $6}')

        bio_rep=$(echo "$bio_raw" | tr -d 'bB')
        tech_rep=$(echo "$tech_raw" | tr -d 'tT')

        cond_lower="$(echo "$cond_raw" | tr '[:upper:]' '[:lower:]')"
        cond_cap="$(echo "${cond_lower:0:1}" | tr '[:lower:]' '[:upper:]')${cond_lower:1}"
        
        target_lower="$(echo "$target_raw" | tr '[:upper:]' '[:lower:]')"
        target_cap="$(echo "${target_lower:0:1}" | tr '[:lower:]' '[:upper:]')${target_lower:1}"

        stype_lower="$(echo "$stype_raw" | tr '[:upper:]' '[:lower:]')"
        if [[ "$stype_lower" == "chip" ]]; then
            sample_type="ChIP"
            target="$target_cap"
            sample_id="${target}_ChIP_${cond_cap}_${time_raw}_rep${bio_rep}"
            input_sample_id="Input_${cond_cap}_${time_raw}_rep1"
        else
            sample_type="Input"
            target="Input"
            sample_id="Input_${cond_cap}_${time_raw}_rep${bio_rep}"
            input_sample_id="NA"
        fi

        # ${time_raw} wird nun in der neuen Spalte ausgegeben
        echo -e "${file_id}\t${sample_id}\t${target}\t${sample_type}\t${input_sample_id}\tWT\t${cond_cap}\t${time_raw}\t${bio_rep}\t${tech_rep}\tfly_head\tSE\t200"
    fi
done < <(sort -u "$tmp_list") >> "$OUTPUT_FILE"

rm -f "$tmp_list"

num_samples=$(($(wc -l < "$OUTPUT_FILE") - 1))
echo "[SUCCESS] '$OUTPUT_FILE' wurde erfolgreich mit ${num_samples} Einträgen generiert!"