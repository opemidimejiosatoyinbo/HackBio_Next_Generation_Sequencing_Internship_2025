```R
#!/bin/rnaseq Rscript
# Script: de_analysis.R
# Description: To perform differential analysis.

# Activate conda environment first

##############################################
# DESeq2 RNA-seq Analysis – Staphylococcus aureus
# Includes enrichment fixes (auto-install, fallback, tryCatch)
##############################################

# Auto-install helper

install_if_missing <- function(pkg, bioc=FALSE) {
if (!requireNamespace(pkg, quietly=TRUE)) {
if (bioc) {
if (!requireNamespace("BiocManager", quietly=TRUE))
install.packages("BiocManager", repos="https://cloud.r-project.org")
BiocManager::install(pkg, ask=FALSE, update=FALSE)
} else {
install.packages(pkg, repos="https://cloud.r-project.org")
}
}
}

# CRAN packages

cran_pkgs <- c("pheatmap","EnhancedVolcano","matrixStats","readr","dplyr","ggplot2")
for (p in cran_pkgs) install_if_missing(p)

# Bioconductor packages

bioc_pkgs <- c("DESeq2","clusterProfiler","KEGGREST","AnnotationDbi")
for (p in bioc_pkgs) install_if_missing(p, bioc=TRUE)

# Load libraries

suppressPackageStartupMessages({
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(matrixStats)
library(readr)
library(dplyr)
library(clusterProfiler)
library(KEGGREST)
library(ggplot2)
library(AnnotationDbi)
})

cat("=== Starting DESeq2 pipeline ===\n")

# Paths

base <- Sys.getenv("HOME")
proj <- file.path(base, "Opemidimeji_2")
counts_file <- file.path(proj, "results/counts/gene_counts_s0.txt")
out_dir <- file.path(proj, "results")
ref_gff3 <- file.path(proj, "ref", "genome_fixed.gff3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Read counts

cat("Loading counts from:", counts_file, "\n")
fc <- read.delim(counts_file, comment.char="#", stringsAsFactors = FALSE, check.names = FALSE)
if (ncol(fc) < 7) stop("featureCounts output looks wrong - not enough columns")
counts <- as.matrix(fc[, 7:ncol(fc)])
rownames(counts) <- fc$Geneid

# Fix sample names

colnames(counts) <- basename(colnames(counts))
colnames(counts) <- sub("\\.sorted\\.bam$", "", colnames(counts))
colnames(counts) <- sub("_Staphylococcus_aureus_RNA-Seq$", "", colnames(counts))
cat("Sample names fixed:\n"); print(colnames(counts))

# Save raw counts

write_csv(as.data.frame(counts), file.path(out_dir, "raw_counts.csv"))
cat("Saved raw_counts.csv\n")

# Metadata

meta_file <- file.path(out_dir, "metadata.csv")
if (!file.exists(meta_file)) {
cat("No metadata.csv found — creating template...\n")
simple_meta <- data.frame(
SampleID = colnames(counts),
Condition = c("chronic","chronic","chronic","chronic","acute","acute","acute")[seq_len(ncol(counts))],
Replicate = rev(seq_len(ncol(counts)))
)
rownames(simple_meta) <- simple_meta$SampleID
write.csv(simple_meta, meta_file, quote=FALSE, row.names=TRUE)
}
metadata <- read.csv(meta_file, row.names = 1, stringsAsFactors = FALSE)
metadata <- metadata[colnames(counts), , drop = FALSE]

# DESeq2

dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ Condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("Condition", "acute", "chronic"))
resOrdered <- res[order(res$padj), ]
res_df <- as.data.frame(resOrdered)
write.csv(res_df, file.path(out_dir, "all_results.csv"))

sig <- subset(res_df, !is.na(padj) & padj < 0.05)
up <- subset(sig, log2FoldChange > 1)
down <- subset(sig, log2FoldChange < -1)
write.csv(up, file.path(out_dir, "upregulated.csv"))
write.csv(down, file.path(out_dir, "downregulated.csv"))

# Variance-stabilizing transform

vsd <- if (nrow(dds) >= 10) vst(dds, blind = FALSE) else rlog(dds, blind = FALSE)

# PCA

png(file.path(out_dir, "PCA.png")); plotPCA(vsd, intgroup = "Condition"); dev.off()

# Volcano

png(file.path(out_dir, "volcano.png"), width=1200, height=900)
EnhancedVolcano(res_df, lab=rownames(res_df), x='log2FoldChange', y='pvalue', pCutoff=0.05, FCcutoff=1)
dev.off()

# Heatmap

top50 <- head(rownames(resOrdered), 50)
if (length(top50) > 1) {
mat <- assay(vsd)[top50, , drop=FALSE]; mat <- mat - rowMeans(mat)
png(file.path(out_dir, "heatmap_top50.png"), width=1000, height=1200)
pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE, annotation_col=metadata)
dev.off()
}

# Save metadata

write.csv(metadata, file.path(out_dir, "metadata_used.csv"))

# Prepare top100 for enrichment

top100 <- sig %>% arrange(padj) %>% head(100)
write.csv(top100, file.path(out_dir, "top100_DE_genes.csv"))

#######################################################

# Enrichment (without org.Sa.eg.db, since not available)

#######################################################
cat("=== Starting enrichment analysis (no org.Sa.eg.db) ===\n")

gene_ids <- rownames(top100)
entrez_ids <- character(0)

# Parse GFF3 to map locus_tag -> GeneID

if (file.exists(ref_gff3)) {
cat("Extracting GeneIDs from GFF3...\n")
gff_lines <- readLines(ref_gff3, warn = FALSE)
locus_to_geneid <- list()
for (ln in gff_lines) {
if (startsWith(ln, "#")) next
f <- strsplit(ln, "\t")[[1]]
if (length(f) < 9) next
attrs <- f[9]
lt <- sub("._locus_tag=([^;]+)._", "\\1", attrs)
if (lt == attrs) lt <- NA
geneid <- NA
m1 <- regmatches(attrs, regexec("GeneID[:=]([0-9]+)", attrs))[[1]]
if (length(m1) >= 2) geneid <- m1[2]
if (!is.na(lt) && !is.na(geneid)) locus_to_geneid[[lt]] <- geneid
}
if (length(locus_to_geneid) > 0) {
mapped <- unlist(locus_to_geneid[intersect(names(locus_to_geneid), gene_ids)])
if (length(mapped) > 0) {
entrez_ids <- unique(mapped)
cat("Mapped via GFF3:", length(entrez_ids), "Entrez IDs found\n")
}
}
}

# GO enrichment (only if IDs available)

if (length(entrez_ids) > 0) {
cat("Running GO enrichment...\n")
ego_bp <- tryCatch(
enrichGO(gene=entrez_ids, OrgDb=NULL, keyType="ENTREZID", ont="BP",
pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE),
error=function(e) NULL
)
if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp))>0) {
write.csv(as.data.frame(ego_bp), file.path(out_dir, "GO_BP_results.csv"))
png(file.path(out_dir, "GO_BP_dotplot.png"), width=1200, height=900)
dotplot(ego_bp, showCategory=20) + ggtitle("GO BP"); dev.off()
} else {
cat("GO enrichment returned no results.\n")
}
} else {
cat("No Entrez IDs → skipping GO enrichment.\n")
}

# KEGG enrichment

if (length(entrez_ids) > 0) {
cat("Running KEGG enrichment with organism='sau'...\n")
kegg_res <- tryCatch(
enrichKEGG(gene=entrez_ids, organism="sau", keyType="ncbi-geneid",
pvalueCutoff=0.05, qvalueCutoff=0.05),
error=function(e) NULL
)
if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res))>0) {
write.csv(as.data.frame(kegg_res), file.path(out_dir, "KEGG_results.csv"))
png(file.path(out_dir, "KEGG_dotplot.png"), width=1200, height=900)
dotplot(kegg_res, showCategory=20) + ggtitle("KEGG"); dev.off()
} else {
cat("KEGG enrichment returned no results.\n")
}
} else {
cat("Skipping KEGG enrichment (no Entrez IDs).\n")
}

cat("=== Pipeline finished successfully ===\n")
```