```bash
#!/bin/bash
set -euo pipefail

# Script: 02_raw_qc.sh

# Description: To run FastQC on raw sequencing data for initial quality assessment.

cd ~/Opemidimeji_2

echo "Starting FastQC on raw FASTQ files"
shopt -s nullglob
for fq in data/raw/\*.fastq.gz; do
base=$(basename "$fq")
echo ">>> Running FastQC on $base"
    fastqc -o results/qc/raw/ "$fq"
echo ">>> Finished FastQC for $base"
done

echo "Running MultiQC summary on raw FastQC reports"
multiqc results/qc/raw/ -o results/qc/raw/

echo "Raw read quality control completed!"
```