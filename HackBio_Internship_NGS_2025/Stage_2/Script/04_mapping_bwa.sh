# Script: 04_mapping_bwa.sh

# Description: To align reads with BWA.

#!/bin/bash
set -euo pipefail

# Directories
BASE_DIR="/home/sararomi/Opemidimeji_2"
DATA_DIR="$BASE_DIR/data"
REF_DIR="$DATA_DIR/ref"
TRIM_DIR="$DATA_DIR/trimmed"
ALIGN_DIR="$BASE_DIR/results/mapping"
LOG_DIR="$BASE_DIR/logs"

mkdir -p "$REF_DIR" "$TRIM_DIR" "$ALIGN_DIR" "$LOG_DIR"

# Correct reference genome URLs (Ensembl Bacteria release 62)
FASTA_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-62/bacteria/fasta/bacteria_113_collection/staphylococcus_aureus_subsp_aureus_usa300_fpr3757_gca_000013465/dna/Staphylococcus_aureus_subsp_aureus_usa300_fpr3757_gca_000013465.ASM1346v1.dna.toplevel.fa.gz"

GFF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-62/bacteria/gff3/bacteria_113_collection/staphylococcus_aureus_subsp_aureus_usa300_fpr3757_gca_000013465/Staphylococcus_aureus_subsp_aureus_usa300_fpr3757_gca_000013465.ASM1346v1.62.gff3.gz"

# File names
GENOME_FA="$REF_DIR/genome.fa"
GENOME_GFF3="$REF_DIR/genome.gff3"

echo ">>> STEP 1: Download reference if missing"

# Download FASTA
if [ ! -s "$GENOME_FA" ]; then
    echo "    - Downloading FASTA..."
    curl -fL "$FASTA_URL" -o "$REF_DIR/genome.fa.gz" || {
        echo "ERROR: FASTA download failed"; exit 1;
    }
    gunzip -c "$REF_DIR/genome.fa.gz" > "$GENOME_FA"
    rm -f "$REF_DIR/genome.fa.gz"
    echo "    - Saved genome FASTA to $GENOME_FA"
else
    echo "    - FASTA already present: $GENOME_FA"
fi

# Download GFF3
if [ ! -s "$GENOME_GFF3" ]; then
    echo "    - Downloading GFF3..."
    curl -fL "$GFF_URL" -o "$REF_DIR/genome.gff3.gz" || {
        echo "ERROR: GFF3 download failed"; exit 1;
    }
    gunzip -c "$REF_DIR/genome.gff3.gz" > "$GENOME_GFF3"
    rm -f "$REF_DIR/genome.gff3.gz"
    echo "    - Saved GFF3 to $GENOME_GFF3"
else
    echo "    - GFF3 already present: $GENOME_GFF3"
fi

echo ">>> STEP 2: Build BWA index"
if [ ! -f "$GENOME_FA.bwt" ]; then
    bwa index "$GENOME_FA"
    echo "    - Index built for $GENOME_FA"
else
    echo "    - BWA index already exists"
fi

echo ">>> STEP 3: Map reads with BWA"
for fq1 in "$TRIM_DIR"/*_1.trim.fastq.gz; do
    fq2=${fq1/_1.trim.fastq.gz/_2.trim.fastq.gz}
    sample=$(basename "$fq1" _1.trim.fastq.gz)
    bam="$ALIGN_DIR/${sample}.sorted.bam"

    if [ ! -s "$bam" ]; then
        echo "    - Mapping $sample"
        bwa mem -t 4 "$GENOME_FA" "$fq1" "$fq2" \
            | samtools view -bS - \
            | samtools sort -o "$bam"
        samtools index "$bam"
    else
        echo "    - BAM already exists: $bam"
    fi
done

echo ">>> STEP 4: Collect mapping stats"
for bam in "$ALIGN_DIR"/*.sorted.bam; do
    samtools flagstat "$bam" > "$LOG_DIR/$(basename "$bam" .sorted.bam).flagstat.txt"
done

echo ">>> Mapping with BWA complete."

for bam in ~/Opemidimeji_2/results/mapping/*.sorted.bam
do
    samtools index $bam
done

# Alternatively before count, redownload NCBI's genome and convert GFF3 to GTF since featureCounts prefers GTF
# To download:
REF_DIR=~/Opemidimeji_2/ref
mkdir -p "$REF_DIR"
cd "$REF_DIR"

BASE="GCF_000013465.1_ASM1346v1"
FASTA_GZ="${BASE}_genomic.fna.gz"
GFF_GZ="${BASE}_genomic.gff.gz"

echo "Downloading genome FASTA..."
curl -fL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/465/${BASE}/${FASTA_GZ}" -o "${FASTA_GZ}"
gunzip -f "${FASTA_GZ}"
mv "${BASE}_genomic.fna" genome.fa

echo "Downloading annotation GFF..."
curl -fL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/465/${BASE}/${GFF_GZ}" -o "${GFF_GZ}"
gunzip -f "${GFF_GZ}"
mv "${BASE}_genomic.gff" genome.gff3

echo "Reference files:"
ls -lh genome.fa genome.gff3

# To convert:
cd ~/Opemidimeji_2/ref

awk 'BEGIN{OFS="\t"} $3=="gene" {
  match($9,/ID=([^;]+)/,a);
  if(a[1]!="") print $1,".","gene",$4,$5,$6,$7,$8,"gene_id \""a[1]"\";"
}' genome.gff3 > genome.gtf