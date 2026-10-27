# ── Setup ─────────────────────────────────────────────────────────────────────
# DESeq2 1.36.0 (R 4.2.2): ecDNA+ (MYCN amp, n = 15) vs ecDNA- (MYCN non-amp, n = 6), design ~ MYCN_amp,
# reference level no_amp, default Wald test / BH FDR. Sample groups come from samples.tsv (audit point 8).
setwd("/path/to/deseq2")                         # folder containing counts_21lines.txt; outputs are written here
SAMPLES_TSV <- "/path/to/repo/samples.tsv"       # repository sample sheet: sample_id, MYCN_amp (amp / no_amp), ...

# ── Load packages ─────────────────────────────────────────────────────────────
library("tidyverse")
library("DESeq2")
library(org.Hs.eg.db)

# ── Load counts ───────────────────────────────────────────────────────────────
# 21-sample matrix from 08_merge_counts.sh step 5 (version-stripped Ensembl IDs).
GE <- read.delim("counts_21lines.txt", row.names = "Geneid", check.names = FALSE)
rownames(GE) <- sub("\\..*$", "", rownames(GE))          # defensive: strip any version suffix
stopifnot(!any(duplicated(rownames(GE))), ncol(GE) == 21)

# ── Load metadata from samples.tsv ────────────────────────────────────────────
samples  <- read.delim(SAMPLES_TSV, stringsAsFactors = FALSE)
stopifnot(all(c("sample_id", "MYCN_amp", "ecDNA_group") %in% colnames(samples)), nrow(samples) == 21)
metadata <- samples[match(colnames(GE), samples$sample_id), c("sample_id", "MYCN_amp", "ecDNA_group", "rna_batch")]
stopifnot(!any(is.na(metadata$sample_id)))                 # every count column must be in samples.tsv
stopifnot(all(metadata$MYCN_amp %in% c("amp", "no_amp")))
stopifnot(all(metadata$ecDNA_group %in% c("ecDNA+", "ecDNA-")))
rownames(metadata) <- metadata$sample_id
metadata$MYCN_amp  <- factor(metadata$MYCN_amp, levels = c("no_amp", "amp"))
# Group definition (audit point 8): the design uses MYCN_amp, the Methods describe ecDNA+ vs ecDNA-. In this 21-line
# cohort the two labels coincide (every MYCN-amplified line carries MYCN on ecDNA), so the contrast is the same.
# Both facts are asserted here rather than assumed.
stopifnot(identical(as.integer(table(metadata$MYCN_amp)), c(6L, 15L)))                   # no_amp = 6, amp = 15
stopifnot(all((metadata$MYCN_amp == "amp") == (metadata$ecDNA_group == "ecDNA+")))     # MYCN_amp == ecDNA_group
message(sprintf("groups: amp/ecDNA+ = %d, no_amp/ecDNA- = %d", sum(metadata$MYCN_amp == "amp"), sum(metadata$MYCN_amp == "no_amp")))
stopifnot(all(colnames(GE) == rownames(metadata)))

# ── Build DESeqDataSet ────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = GE,
                              colData   = metadata,
                              design    = ~MYCN_amp)

# Reference level: no_amp = ecDNA− (denominator); log2FC > 0 means higher in ecDNA+ (amp)
dds$MYCN_amp <- relevel(dds$MYCN_amp, ref = "no_amp")

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

# ── Session record ────────────────────────────────────────────────────────────
writeLines(capture.output(sessionInfo()), "deseq2_sessionInfo.txt")
