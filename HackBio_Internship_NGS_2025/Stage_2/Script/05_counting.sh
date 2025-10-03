# Script: 05_counting.sh

# Description: To use featureCounts to generate a count matrix for DESeq2.

featureCount to generate a count matrix for DESeq2.

#!/bin/bash

# Using GFF3 directly (no gffread)
# Set working directories
WORKDIR=/home/sararomi/Opemidimeji_2
REF=$WORKDIR/ref
MAP=$WORKDIR/results/mapping
COUNT=$WORKDIR/results/counts

mkdir -p $COUNT

echo ">>> Step 1: Run featureCounts directly with GFF3"

# Unstranded (-s 0)
featureCounts -T 8 -p -F GFF3 -t gene -g locus_tag \
  -a $REF/genome_fixed.gff3 \
  -o $COUNT/gene_counts_s0.txt \
  $MAP/*.sorted.bam

# Forward stranded (-s 1)
featureCounts -T 8 -p -F GFF3 -t gene -g locus_tag -s 1 \
  -a $REF/genome_fixed.gff3 \
  -o $COUNT/gene_counts_s1.txt \
  $MAP/*.sorted.bam

# Reverse stranded (-s 2)
featureCounts -T 8 -p -F GFF3 -t gene -g locus_tag -s 2 \
  -a $REF/genome_fixed.gff3 \
  -o $COUNT/gene_counts_s2.txt \
  $MAP/*.sorted.bam

echo ">>> Step 2: Summaries"
ls -lh $COUNT/*.summary

echo ">>> Step 3: Compare Assigned reads for strandedness choice"
echo -e "Strandedness\tAssigned_reads"
awk '/Assigned/ {print "s0\t"$2}' $COUNT/gene_counts_s0.txt.summary
awk '/Assigned/ {print "s1\t"$2}' $COUNT/gene_counts_s1.txt.summary
awk '/Assigned/ {print "s2\t"$2}' $COUNT/gene_counts_s2.txt.summary

echo "featureCounts completed!"