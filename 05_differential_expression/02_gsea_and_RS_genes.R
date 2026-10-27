library(tidyverse)
library(fgsea)
library(msigdbr)
library(ggplot2)

# dplyr verbs are written as dplyr::<verb>() where a Bioconductor package attached
# earlier in the same session (DESeq2 -> S4Vectors, org.Hs.eg.db -> AnnotationDbi)
# would otherwise mask them: `rename` and `select`. Running this script in a fresh
# R session is still recommended.

# ── Paths ──────────────────────────────────────────────────────────────────────
deseq2_path <- '/path/to/deseq2/deseq2_results_ecDNA.csv'   # output of 01_deseq2_ecDNA.R
out_dir     <- '/path/to/gsea'

# ── Versions ──────────────────────────────────────────────────────────────────
# Gene-set membership changes between MSigDB releases, so Hallmark results are tied
# to the msigdbr version below. Original analysis: msigdbr 26.1.0 (MSigDB 2026.1.Hs),
# fgsea 1.22.0, R 4.2.2 -- see README "Software versions" and the per-run log written here.
stopifnot(packageVersion("msigdbr") >= "26.1.0")
writeLines(c(
  paste("R      ", R.version.string),
  paste("fgsea  ", packageVersion("fgsea")),
  paste("msigdbr", packageVersion("msigdbr")),
  paste("dplyr  ", packageVersion("dplyr")),
  paste("ggplot2", packageVersion("ggplot2")),
  paste("run    ", format(Sys.time(), "%Y-%m-%d %H:%M"))
), file.path(out_dir, "gsea_session_versions.txt"))

# ── Figure text scaling ───────────────────────────────────────────────────────
# Change only this value to enlarge or reduce all text in both figures.
FONT_SCALE <- 1.9

# ── Step 1: Build ranked gene list ────────────────────────────────────────────
# Column 1 of the DESeq2 CSV is the unnamed row-name column (Ensembl ID).
df <- read.csv(deseq2_path) %>%
  dplyr::rename(ensembl = 1) %>%
  filter(!is.na(pvalue), !is.na(log2FoldChange), !is.na(Symbol)) %>%
  filter(!duplicated(Symbol)) %>%
  filter(pvalue < 1) %>%
  mutate(
    rank = sign(log2FoldChange) * -log10(pvalue)
    + log2FoldChange * 1e-5
  ) %>%
  arrange(desc(rank))

ranked <- setNames(df$rank, df$Symbol)
cat("Genes in ranked list:", length(ranked), "\n")

# ── Step 2: Get gene sets ──────────────────────────────────────────────────────
# `collection =` replaces the deprecated `category =` argument (msigdbr >= 10).
hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

cat("Hallmark gene sets:", length(hallmark), "\n")

# ── Step 3: Run fgsea on Hallmark ─────────────────────────────────────────────
# fgseaMultilevel is stochastic; the seed is set inside the function so every call
# is reproducible on its own (seed 42, as in the original run).
run_fgsea <- function(pathways, stats, label, seed = 42) {
  set.seed(seed)
  fgseaMultilevel(
    pathways = pathways,
    stats    = stats,
    minSize  = 15,
    maxSize  = 500,
    eps      = 0
  ) %>%
    arrange(desc(NES)) %>%
    mutate(collection = label)
}

results_hallmark <- run_fgsea(hallmark, ranked, "Hallmark")

write.csv(results_hallmark %>% dplyr::select(-leadingEdge),
          file.path(out_dir, "hallmark_GSEA_results.csv"), row.names = FALSE)

cat("Hallmark:", nrow(results_hallmark), "pathways tested;",
    sum(results_hallmark$padj < 0.05), "with FDR < 0.05\n")

# ── Step 4: Hallmark dot plot — top 20 significant non-MYC sets ───────────────
# MYC_TARGETS_V1/V2 are excluded from the plot (expected given MYCN amplification);
# they remain in hallmark_GSEA_results.csv.
results_hallmark <- results_hallmark %>%
  mutate(is_myc = grepl("MYC", pathway, ignore.case = TRUE))

non_myc <- filter(results_hallmark, !is_myc)

plot_df_hallmark <- non_myc %>%
  filter(padj < 0.05) %>%                        # keep only significant
  slice_max(abs(NES), n = 20) %>%                 # top 20 by absolute NES
  arrange(NES) %>%
  mutate(
    pathway_clean = str_replace(pathway, "HALLMARK_", "") %>%
      str_replace_all("_", " "),
    direction     = ifelse(NES > 0, "Enriched in ecDNA+", "Enriched in ecDNA-"),
    pathway_clean = factor(pathway_clean, levels = pathway_clean)
  )

p_hallmark <- ggplot(plot_df_hallmark, aes(x = NES, y = pathway_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.5) +
  geom_point(aes(color = direction, size = -log10(pmax(padj, 1e-4))), alpha = 0.85) +
  scale_color_manual(values = c("Enriched in ecDNA+" = "#C0392B",
                                "Enriched in ecDNA-" = "#2471A3")) +
  guides(
    color = guide_legend(
      order = 1, nrow = 1,
      override.aes = list(size = 5)
    ),
    size = guide_legend(
      order = 2, nrow = 2, byrow = TRUE,
      title.position = "top"
    )
  ) +
  scale_size_continuous(name = "-log10(FDR)", range = c(2, 8)) +
  labs(
    title    = NULL,
    subtitle = NULL,
    x = "Normalized Enrichment Score (NES)", y = NULL, color = NULL
  ) +
  theme_minimal(base_size = 10 * FONT_SCALE) +
  theme(
    panel.grid      = element_blank(),
    axis.line.x     = element_line(color = "grey80"),
    legend.position = "top",
    legend.box      = "vertical",   # stack color and size guides to prevent clipping
    legend.direction = "horizontal",
    legend.box.just = "center",
    legend.spacing.y = grid::unit(1, "pt"),
    legend.key.spacing.x = grid::unit(3, "pt"),
    legend.margin   = margin(2, 2, 2, 2),
    plot.margin     = margin(5, 20, 20, 5),
    legend.text     = element_text(size = 14),
    legend.title    = element_text(size = 14),
    axis.title.x = element_text(size = 8 * FONT_SCALE),
    axis.text.y     = element_text(size = 8 * FONT_SCALE),
    axis.text.x     = element_text(size = 8 * FONT_SCALE)
  )

ggsave(file.path(out_dir, "GSEA_hallmark_dotplot.pdf"),
       p_hallmark, width = 10, height = 8, dpi = 900)
cat("Hallmark dot plot saved.\n")

# ── Step 5: RS gene bar plot ──────────────────────────────────────────────────
# Curated replication-stress (RS) gene panel used for Fig 2E (n = length(rs_genes)).
# This is the only RS list in the repository; 03_export_CHEK1_FPKM.R reads the same
# vector from rs_genes.txt so the two scripts cannot drift apart.
rs_genes <- c(
  # ATR-CHK1 axis
  'ATR', 'ATRIP', 'TOPBP1', 'CHEK1', 'CLSPN', 'TIMELESS', 'TIPIN', 'RAD17', 'WEE1',
  # Parallel checkpoint
  'CHEK2',
  # RPA complex
  'RPA1', 'RPA2', 'RPA3',
  # Nucleotide pool
  'RRM1', 'RRM2',
  # Replisome
  'MCM2', 'MCM7', 'CDC45', 'PCNA', 'RFC1', 'POLD1', 'POLE', 'CDC6',
  # Fork processing / nucleases
  'FEN1', 'MRE11',
  # HR repair / fork protection
  'RAD51', 'BRCA1', 'BRCA2',
  # Fanconi anemia
  'FANCD2'
)
writeLines(rs_genes, file.path(out_dir, "rs_genes.txt"))

rs_df <- df %>%
  filter(Symbol %in% rs_genes) %>%
  arrange(log2FoldChange) %>%
  mutate(
    Symbol    = factor(Symbol, levels = Symbol),
    # Bar colour: exploratory threshold FDR < 0.2 (stated as such in the Fig 2E legend).
    # Asterisks: conventional FDR < 0.05 / 0.01 / 0.001.
    sig_color = case_when(
      padj < 0.2 & log2FoldChange > 0 ~ "Higher in ecDNA+ (FDR<0.2)",
      padj < 0.2 & log2FoldChange < 0 ~ "Higher in ecDNA- (FDR<0.2)",
      TRUE                            ~ "Not significant"
    ),
    sig_label = case_when(
      padj < 0.001 ~ "***", padj < 0.01 ~ "**",
      padj < 0.05  ~ "*",  TRUE ~ ""
    )
  )

mean_fc <- mean(rs_df$log2FoldChange)
cat(sprintf("Mean log2FC of RS genes: %.3f\n", mean_fc))
cat(sprintf("RS genes found in DESeq2 results: %d / %d\n", nrow(rs_df), length(rs_genes)))
cat(sprintf("RS genes with FDR < 0.05: %s\n",
            paste(as.character(rs_df$Symbol[rs_df$padj < 0.05]), collapse = ", ")))

# Report any RS genes missing from the DESeq2 results
missing <- setdiff(rs_genes, rs_df$Symbol)
if (length(missing) > 0) {
  cat("RS genes not found in DESeq2 results:", paste(missing, collapse = ", "), "\n")
}

p_rs <- ggplot(rs_df, aes(x = log2FoldChange, y = Symbol, fill = sig_color)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.5) +
  geom_col(color = "white", linewidth = 0.3, width = 0.65) +
  geom_text(aes(label = sig_label,
                x     = log2FoldChange + ifelse(log2FoldChange >= 0, 0.08, -0.08),
                hjust = ifelse(log2FoldChange >= 0, 0, 1)),
            size = 3.5 * FONT_SCALE, fontface = "bold", color = "#2C3E50") +
  scale_fill_manual(values = c(
    "Higher in ecDNA+ (FDR<0.2)" = "#C0392B",
    "Higher in ecDNA- (FDR<0.2)" = "#2471A3",
    "Not significant"            = "#D5D8DC"
  )) +
  labs(
    title    = NULL,
    subtitle = NULL,
    x = expression(log[2]~"fold change  (ecDNA+ / ecDNA-)"),
    y = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 10 * FONT_SCALE) +
  theme(
    panel.grid      = element_blank(),
    axis.line.x     = element_line(color = "grey80"),
    legend.position = "top"
  )

ggsave(file.path(out_dir, "RS_genes_log2FC.pdf"),
       p_rs, width = 7, height = 10, dpi = 900)
cat("RS gene plot saved.\n")

# ── Step 6: Summary ───────────────────────────────────────────────────────────
cat("\n================================================\n")
cat("OUTPUT FILES\n")
cat("================================================\n")
cat("hallmark_GSEA_results.csv   -- all 50 Hallmark sets, fgsea statistics\n")
cat("GSEA_hallmark_dotplot.pdf   -- top 20 significant non-MYC sets by |NES| (Fig 2D)\n")
cat(sprintf("RS_genes_log2FC.pdf         -- RS gene bar plot (n=%d), colored = FDR<0.2 (Fig 2E)\n",
            length(rs_genes)))
cat("rs_genes.txt                -- the RS gene panel, one symbol per line\n")
cat("gsea_session_versions.txt   -- R / fgsea / msigdbr / dplyr / ggplot2 versions for this run\n")
