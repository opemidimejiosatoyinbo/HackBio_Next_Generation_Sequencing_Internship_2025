#!/bin/bash
   # Script: downsampling.sh
   # Description: To choose 50 samples randomly from the raw_data

   # Set directories
   PARENT_DIR="Opemidimeji/data/raw_data"
   RAW_DATA_50="Opemidimeji/data/raw_data_50"

   # Create a destination directory
   mkdir -p "$RAW_DATA_50"

   # Get a list of all the R1 files and shuffle them randomly, then pick the first 50
   find "$PARENT_DIR" -maxdepth 1 -name "*_1.fastq.gz" | shuf | head -50 | while read r1_file; do
       mv "$r1_file" "$RAW_DATA_50/"

       # Find and move the corresponding R2 file
       base_name=$(basename "$r1_file" _1.fastq.gz)
       r2_file="$PARENT_DIR/${base_name}_2.fastq.gz"
       mv "$r2_file" "$RAW_DATA_50/"
   done