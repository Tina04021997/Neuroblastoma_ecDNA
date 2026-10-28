#!/usr/bin/env Rscript
# Sensitivity analysis (not used for any figure): does adding the RNA-seq sequencing batch as a covariate change the
# ecDNA+ vs ecDNA- DESeq2 signature? The published model is ~MYCN_amp (01_deseq2_ecDNA.R); the two batches are not
# confounded with the groups (batch 1: 6 amp / 4 no_amp; batch 2: 9 amp / 2 no_amp), so ~rna_batch + MYCN_amp is estimable.
# Usage: Rscript 04_batch_sensitivity.R <counts_21lines.txt> <samples.tsv> <deseq2_results_ecDNA.csv> [out.tsv]
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: Rscript 04_batch_sensitivity.R <counts_21lines.txt> <samples.tsv> <deseq2_results_ecDNA.csv> [out.tsv]")
GE  <- read.delim(args[1], row.names = "Geneid", check.names = FALSE)
s   <- read.delim(args[2]); s <- s[match(colnames(GE), s$sample_id), ]
stopifnot(!any(is.na(s$sample_id)), nrow(s) == 21)
pub <- read.csv(args[3], row.names = 1)                       # published universe: symbol-mapped, duplicate symbols dropped
md  <- data.frame(row.names = s$sample_id, MYCN_amp = factor(s$MYCN_amp, levels = c("no_amp", "amp")), batch = factor(s$rna_batch))
print(table(MYCN_amp = md$MYCN_amp, batch = md$batch))
run <- function(design) {
  dds <- DESeq(DESeqDataSetFromMatrix(GE, md, design), quiet = TRUE)
  as.data.frame(results(dds, name = "MYCN_amp_amp_vs_no_amp"))[rownames(pub), ]
}
r0 <- run(~MYCN_amp); r1 <- run(~batch + MYCN_amp)
stopifnot(isTRUE(all.equal(r0$padj, pub$padj)))                 # the unadjusted model must reproduce the published export
sig <- function(r) { r <- r[!is.na(r$padj), ]; list(up = rownames(r)[r$padj < 0.05 & r$log2FoldChange > 1], down = rownames(r)[r$padj < 0.05 & r$log2FoldChange < -1]) }
a <- sig(r0); b <- sig(r1)
common <- rownames(r0)[!is.na(r0$padj) & !is.na(r1$padj)]
cat(sprintf("~MYCN_amp          : %d up / %d down\n~batch + MYCN_amp  : %d up / %d down\n", length(a$up), length(a$down), length(b$up), length(b$down)))
cat(sprintf("published DE genes still DE after batch adjustment: %d of %d (%.0f%%); log2FC Pearson r = %.3f (%d genes)\n",
            length(intersect(a$up, b$up)) + length(intersect(a$down, b$down)), length(a$up) + length(a$down),
            100 * (length(intersect(a$up, b$up)) + length(intersect(a$down, b$down))) / (length(a$up) + length(a$down)),
            cor(r0[common, "log2FoldChange"], r1[common, "log2FoldChange"]), length(common)))
for (g in c("MYCN", "CHEK1", "MCM7", "FANCD2", "RRM1", "CDC45", "CDC6")) { i <- which(pub$Symbol == g)
  cat(sprintf("  %-7s log2FC %5.2f -> %5.2f   padj %.3f -> %.3f\n", g, r0[i, "log2FoldChange"], r1[i, "log2FoldChange"], r0[i, "padj"], r1[i, "padj"])) }
if (length(args) >= 4) {
  out <- data.frame(Ensembl = rownames(pub), Symbol = pub$Symbol, log2FC_unadjusted = r0$log2FoldChange, padj_unadjusted = r0$padj,
                    log2FC_batch_adjusted = r1$log2FoldChange, padj_batch_adjusted = r1$padj)
  write.table(out, args[4], sep = "\t", quote = FALSE, row.names = FALSE); cat("written:", args[4], "\n")
}
