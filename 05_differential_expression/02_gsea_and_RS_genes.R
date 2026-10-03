library(tidyverse)
library(fgsea)
library(msigdbr)
library(ggplot2)

# ── Paths ──────────────────────────────────────────────────────────────────────
deseq2_path <- '/path/to/deseq2/deseq2_results_ecDNA.csv'   # output of 01_deseq2_ecDNA.R
out_dir     <- '/path/to/gsea'

# ── Figure text scaling ───────────────────────────────────────────────────────
# Change only this value to enlarge or reduce all text in both figures.
FONT_SCALE <- 1.9

# ── Step 1: Build ranked gene list ────────────────────────────────────────────
df <- read.csv(deseq2_path) %>%
  rename(ensembl = 1) %>%
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
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

cat("Hallmark gene sets:", length(hallmark), "\n")

# ── Step 3: Run fgsea on Hallmark ─────────────────────────────────────────────
set.seed(42)

run_fgsea <- function(pathways, stats, label) {
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

write.csv(results_hallmark %>% select(-leadingEdge),
          file.path(out_dir, "hallmark_GSEA_results.csv"), row.names = FALSE)

cat("Hallmark:", nrow(results_hallmark), "pathways tested\n")

# ── Step 4: Hallmark dot plot — top 20 significant, no red labels ──────────────
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

# ── Step 5: RS gene bar plot — expanded gene list ─────────────────────────────
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

rs_df <- df %>%
  filter(Symbol %in% rs_genes) %>%
  arrange(log2FoldChange) %>%
  mutate(
    Symbol    = factor(Symbol, levels = Symbol),
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
cat("hallmark_GSEA_results.csv\n")
cat("GSEA_hallmark_dotplot.pdf    -- Hallmark top 20 significant, no red labels\n")
cat("RS_genes_log2FC.pdf         -- RS gene bar plot (n=30), colored = FDR<0.2\n")
