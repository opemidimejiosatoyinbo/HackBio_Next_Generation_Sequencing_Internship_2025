```bash
#!/bin/bash
set -euo pipefail
shopt -s nullglob

# Script: 03_trim_and_qc.sh

# Description: Perform trimming with fastp followed by quality control (FastQC + MultiQC) for single-end (transcriptome) FASTQ

cd ~/Opemidimeji_2

echo "Starting trimming and QC with fastp"

# collect R1

r1_files=(data/raw/_\_1.fastq.gz)
if [ ${#r1_files[@]} -eq 0 ]; then
echo "!!! ERROR: No files found matching _\_1.fastq.gz in data/raw"
ls -lh data/raw | head -n 20
exit 1
fi

for fq1 in "${r1_files[@]}"; do
    base=$(basename "$fq1" _1.fastq.gz)
    fq2="data/raw/${base}\_2.fastq.gz"

    echo ">>> Processing sample: $base"
    if [ ! -f "$fq2" ]; then
        echo "!!! ERROR: Paired file not found for $fq1 (expected $fq2)"
        continue
    fi

    fastp \
      -i "$fq1" -I "$fq2" \
      -o "data/trimmed/${base}_1.trim.fastq.gz" \
      -O "data/trimmed/${base}_2.trim.fastq.gz" \
      --detect_adapter_for_pe \
      --cut_front --cut_tail --cut_mean_quality 20 \
      --length_required 50 \
      --html "results/qc/trimmed/${base}_fastp.html" \
      --json "results/qc/trimmed/${base}_fastp.json" \
      --thread 8

    if [ $? -eq 0 ]; then
        echo ">>> Finished trimming $base"
        echo "    Output: data/trimmed/${base}_1.trim.fastq.gz, data/trimmed/${base}_2.trim.fastq.gz"
        echo "    QC reports: results/qc/trimmed/${base}_fastp.html, results/qc/trimmed/${base}_fastp.json"
    else
        echo "!!! ERROR: fastp failed for $base"
    fi

done

# Optionally run MultiQC on trimmed fastp HTMLs

if compgen -G "results/qc/trimmed/\*.html" > /dev/null; then
multiqc results/qc/trimmed -o results/qc/trimmed
fi

echo "Trimming completed for all samples!"
```