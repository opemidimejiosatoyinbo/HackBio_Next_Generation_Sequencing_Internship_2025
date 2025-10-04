### 1. Create Directory

```bash
#!/bin/bash

# Create project directory structure

set -euo pipefail

BASE=~/Opemidimeji_2

mkdir -p "$BASE"/{scripts,data/{raw,trimmed,aln,ref},results/{qc/raw,qc/trimmed,counts,mapping,plots}}
echo "Project directory structure created under $BASE"
```

### 2. Download

Navigate to 'Opemidimeji_2/data/raw', Use SRA-Explorer to download the samples (4 chronic PJI and 3 acute PJI) from PRJNA867318

```bash
#!/bin/bash
set -euo pipefail

# Script: 01_download.sh

# Description: To download samples (4 chronic PJI and 4 acute PJI) from PRJNA867318

# Download raw FASTQ files (paired-end) from SRA into Opemidimeji_2/data/raw

cd ~/Opemidimeji_2/data/raw

echo ">>> STEP: Downloading FASTQ files to $(pwd)"

# Chronic PJI samples

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/077/SRR20959677/SRR20959677_1.fastq.gz -o SRR20959677_GSM6435907_C-PJI03_Staphylococcus_aureus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/077/SRR20959677/SRR20959677_2.fastq.gz -o SRR20959677_GSM6435907_C-PJI03_Staphylococcus_aureus_RNA-Seq_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/079/SRR20959679/SRR20959679_1.fastq.gz -o SRR20959679_GSM6435905_C-PJI01_Staphylococcus_aureus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/079/SRR20959679/SRR20959679_2.fastq.gz -o SRR20959679_GSM6435905_C-PJI01_Staphylococcus_aureus_RNA-Seq_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/078/SRR20959678/SRR20959678_1.fastq.gz -o SRR20959678_GSM6435906_C-PJI02_Staphylococcus_aureus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/078/SRR20959678/SRR20959678_2.fastq.gz -o SRR20959678_GSM6435906_C-PJI02_Staphylococcus_aureus_RNA-Seq_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/076/SRR20959676/SRR20959676_1.fastq.gz -o SRR20959676_GSM6435908_C-PJI04_Staphylococcus_aureus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/076/SRR20959676/SRR20959676_2.fastq.gz -o SRR20959676_GSM6435908_C-PJI04_Staphylococcus_aureus_RNA-Seq_2.fastq.gz

# Acute PJI samples

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/082/SRR20959682/SRR20959682_1.fastq.gz -o SRR20959682_GSM6435902_A-PJI02_Staphylococcus_aureus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/082/SRR20959682/SRR20959682_2.fastq.gz -o SRR20959682_GSM6435902_A-PJI02_Staphylococcus_aureus_RNA-Seq_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/081/SRR20959681/SRR20959681_1.fastq.gz -o SRR20959681_GSM6435903_A-PJI03_Staphylococcus_aureus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/081/SRR20959681/SRR20959681_2.fastq.gz -o SRR20959681_GSM6435903_A-PJI03_Staphylococcus_aureus_RNA-Seq_2.fastq.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/080/SRR20959680/SRR20959680_1.fastq.gz -o SRR20959680_GSM6435904_A-PJI04_Staphylococcus_aureus_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/080/SRR20959680/SRR20959680_2.fastq.gz -o SRR20959680_GSM6435904_A-PJI04_Staphylococcus_aureus_RNA-Seq_2.fastq.gz

echo "All samples downloaded!"
```