#!/usr/bin/env Rscript
# ============================================================
# RS Gene Expression Heatmap
# ============================================================
# Dependencies:
#   install.packages(c("tidyverse","pheatmap","RColorBrewer"))
# Expects: ${sample}.gene.abundance.txt files in working directory
# ============================================================

setwd("/path/to/stringtie")   # folder containing <sample>.gene.abundance.txt (StringTie); outputs are written here

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

# ── Sample metadata ──────────────────────────────────────────
ecdna <- c(
  IMR32 = 1, Kelly = 1, NGP = 2, LAN6 = 0, NBLS = 0,
  TR14 = 5, CHP212 = 1, CHP134 = 1, SKNAS = 0, SKNSH = 0,
  KANR = 1, KAN = 1, KCNR = 1, KPNYN = 1, LAN1 = 2,
  LANS = 1, NB1643 = 1, NB69 = 0, NB10 = 1, SKNBE2 = 1,
  SKNrI = 0
)

samples     <- names(ecdna)
ecdna_group <- ifelse(ecdna > 0, "ecDNA+", "ecDNA-")

# ── RS genes of interest ─────────────────────────────────────
rs_genes <- c("FEN1","ATR","CHEK1","RRM2","RPA1",
              "RPA2","RAD51","BRCA1","CDC6","CLSPN")

# ── Load all sample files ─────────────────────────────────────
message("Loading gene abundance files...")

load_sample <- function(sample) {
  fname <- paste0(sample, ".gene.abundance.txt")
  if (!file.exists(fname)) {
    warning("File not found: ", fname)
    return(NULL)
  }
  df <- read.table(fname, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
  df$Sample <- sample
  df
}

all_data <- map(samples, load_sample) |>
  compact() |>
  bind_rows()

loaded_samples <- intersect(samples, unique(all_data$Sample))
message(sprintf("Loaded %d / %d samples.", length(loaded_samples), length(samples)))

# ── Filter to RS genes and transform ─────────────────────────
rs_data <- all_data |>
  filter(Gene.Name %in% rs_genes) |>
  select(Sample, Gene.Name, FPKM) |>
  mutate(
    log2FPKM    = log2(FPKM + 1),
    ecDNA_group = ecdna_group[Sample]
  )

# Guard: warn about missing genes
missing <- setdiff(rs_genes, unique(rs_data$Gene.Name))
if (length(missing)) warning("Genes not found in data: ", paste(missing, collapse = ", "))

# Guard: catch sample name mismatches before silent NAs
na_samples <- rs_data |> filter(is.na(ecDNA_group)) |> pull(Sample) |> unique()
if (length(na_samples)) {
  stop("Sample names don't match ecdna dictionary - check capitalisation: ",
       paste(na_samples, collapse = ", "))
}

# ── Build gene × sample matrix ───────────────────────────────
mat <- rs_data |>
  select(Sample, Gene.Name, log2FPKM) |>
  pivot_wider(names_from = Sample, values_from = log2FPKM,
              values_fill = NA) |>
  column_to_rownames("Gene.Name") |>
  as.matrix()

# ── Column annotation (ecDNA status) ─────────────────────────
col_anno <- data.frame(
  ecDNA     = ecdna_group[colnames(mat)],
  row.names = colnames(mat)
)

anno_colors <- list(
  ecDNA = c("ecDNA+" = "#E64B35", "ecDNA-" = "#4DBBD5")
)

# Sort columns: ecDNA+ first, then ecDNA-
col_order       <- order(col_anno$ecDNA, decreasing = TRUE)
mat_sorted      <- mat[, col_order]
col_anno_sorted <- as.data.frame(col_anno[col_order, , drop = FALSE])
colnames(col_anno_sorted) <- "ecDNA"  # preserve column name after row subset

# ── Draw and save heatmap ─────────────────────────────────────
message("Generating heatmap...")

png("RS_heatmap.png", width = 10, height = 6, units = "in", res = 300)
pheatmap(
  mat_sorted,
  annotation_col    = col_anno_sorted,
  annotation_colors = anno_colors,
  cluster_cols      = FALSE,          # columns pre-sorted by ecDNA status
  cluster_rows      = TRUE,           # cluster genes by expression pattern
  scale             = "row",          # z-score each gene across all samples
  color             = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
  border_color      = NA,
  fontsize          = 9,
  main              = "RS Gene Expression (row z-score of log2[FPKM+1])",
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  angle_col         = 45
)
dev.off()

message("Saved: Fig4A_heatmap.png")
message("\n✓ Done.")

# ── Print CHEK1 FPKM for all samples ─────────────────────────
chek1_fpkm <- rs_data |>
  filter(Gene.Name == "CHEK1") |>
  select(Sample, ecDNA_group, FPKM) |>
  arrange(desc(ecDNA_group), Sample)   # ecDNA+ first, then ecDNA-

print(chek1_fpkm, n = Inf)

# Optional: also write to a file
write.csv(chek1_fpkm, "CHEK1_FPKM.csv", row.names = FALSE)