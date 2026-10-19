# ── Setup ─────────────────────────────────────────────────────────────────────
setwd("/path/to/deseq2")   # folder containing counts_21lines.txt and heatmap_condition.csv; outputs are written here

# ── Load packages ─────────────────────────────────────────────────────────────
library("tidyverse")
library("DESeq2")
library(org.Hs.eg.db)

# ── Load counts ───────────────────────────────────────────────────────────────
# 21-sample matrix from 08_merge_counts.sh step 5 (version-stripped Ensembl IDs).
GE <- read.delim("counts_21lines.txt", row.names = "Geneid", check.names = FALSE)
rownames(GE) <- sub("\\..*$", "", rownames(GE))          # defensive: strip any version suffix
stopifnot(!any(duplicated(rownames(GE))), ncol(GE) == 21)

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
res_df         <- as.data.frame(res)
res_df$Ensembl <- rownames(res_df)
res_df$Symbol  <- mapIds(org.Hs.eg.db, res_df$Ensembl,
                         keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")

n_mapped <- sum(!is.na(res_df$Symbol))
message(sprintf("Symbol mapping: %d / %d genes", n_mapped, nrow(res_df)))
stopifnot(n_mapped / nrow(res_df) > 0.5)                  # would be ~0 with versioned IDs

 # Duplicate symbols: keep the first occurrence in file (GTF) order — deterministic for a fixed annotation and identical to the original run. Report how many were dropped.
res_df <- res_df[!is.na(res_df$Symbol), ]
n_dup  <- sum(duplicated(res_df$Symbol))
res_df <- res_df[!duplicated(res_df$Symbol), ]
message(sprintf("Dropped %d duplicate-symbol rows", n_dup))

write.csv(res_df, file = "deseq2_results_ecDNA.csv")   # Ensembl ID stays as row name (GSEA script reads column 1)
