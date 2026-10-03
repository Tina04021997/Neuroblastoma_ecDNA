# ── Setup ─────────────────────────────────────────────────────────────────────
setwd("/path/to/deseq2")   # folder containing merged_counts.txt and heatmap_condition.csv; outputs are written here

# ── Load packages ─────────────────────────────────────────────────────────────
library("tidyverse")
library("DESeq2")
library(org.Hs.eg.db)

# ── Load counts ───────────────────────────────────────────────────────────────
GE <- read.delim("merged_counts.txt", row.names="Geneid")
dim(GE)

# ── Load metadata ─────────────────────────────────────────────────────────────
metadata <- read_csv("./heatmap_condition.csv")
metadata <- as.data.frame(metadata)

# Sanity check — column order must match metadata order
stopifnot(all(colnames(GE) == metadata$Sample))

# ── Build DESeqDataSet ────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = GE,
                              colData   = metadata,
                              design    = ~MYCN_amp)

# Set reference level: no_amp = ecDNA− (denominator)
# log2FC > 0 means higher in ecDNA+ (amp)
dds$MYCN_amp <- relevel(dds$MYCN_amp, ref="no_amp")

# ── Run DESeq2 ────────────────────────────────────────────────────────────────
dds <- DESeq(dds)
res <- results(dds)
summary(res)

# ── Add gene symbols and export ───────────────────────────────────────────────
res_df        <- as.data.frame(res)
res_df$Symbol <- mapIds(org.Hs.eg.db,
                        rownames(res_df),
                        keytype = "ENSEMBL",
                        column  = "SYMBOL")

# Keep only genes with a symbol, remove duplicates
res_df <- res_df[!is.na(res_df$Symbol), ]
res_df <- res_df[!duplicated(res_df$Symbol), ]

write.csv(res_df, file="deseq2_results_ecDNA.csv")

print(paste("Total genes exported:", nrow(res_df)))
print(paste("Significant (FDR < 0.05, |log2FC| > 1):",
            sum(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05, na.rm=TRUE)))