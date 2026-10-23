#!/usr/bin/env Rscript
# ============================================================
# Export StringTie FPKM tables used by the figures
# ============================================================
# Reads <sample>.gene.abundance.txt (StringTie -A output, 06_stringtie.sh) for the
# 21 lines listed in samples.tsv and writes:
#   MYCN_FPKM.txt        Sample / MYCN / MYCN_amp          -> 06_figures/fig2B_MYCN_FPKM_boxplot.py
#   MYCN_CHEK1_FPKM.csv  Sample / ecDNA_group / MYCN / CHEK1 -> MYCN_FPKM and CHEK1_FPKM columns of
#                                                             data/master.csv (Fig 1A oncoprint)
#   CHEK1_FPKM.csv       Sample / ecDNA_group / FPKM        (kept for backward compatibility)
# and, optionally, an exploratory RS-gene heatmap (RS_heatmap.png; not a manuscript figure).
#
# Sample groups come from samples.tsv (no hard-coded sample dictionary here), and the RS
# panel is read from rs_genes.txt written by 02_gsea_and_RS_genes.R, so the repository
# carries exactly one definition of each.
#
# Dependencies: tidyverse, pheatmap, RColorBrewer
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

# ── Paths: edit before running ───────────────────────────────────────────────
STRINGTIE_DIR <- "/Users/tingyang/Desktop/LAB/Curebound/Pilot_Second/RS"          # folder containing <sample>.gene.abundance.txt
SAMPLES_TSV   <- "/Users/tingyang/Desktop/LAB/Curebound/Pilot_Second/FPKM_2026_rerun/samples.tsv"
RS_GENES_TXT  <- "/Users/tingyang/Desktop/LAB/Curebound/Pilot_Second/GSEA/2026_rerun/rs_genes.txt"  # written by 02_gsea_and_RS_genes.R
OUT_DIR       <- STRINGTIE_DIR
DRAW_HEATMAP  <- TRUE
# ─────────────────────────────────────────────────────────────────────────────

# ── Sample metadata ──────────────────────────────────────────
samples <- read_tsv(SAMPLES_TSV, show_col_types = FALSE)
stopifnot(all(c("sample_id", "MYCN_amp", "ecDNA_group") %in% names(samples)),
          nrow(samples) == 21)
ecdna_group <- setNames(samples$ecDNA_group, samples$sample_id)

# ── Load all sample files ─────────────────────────────────────
message("Loading gene abundance files...")

load_sample <- function(sample) {
  fname <- file.path(STRINGTIE_DIR, paste0(sample, ".gene.abundance.txt"))
  if (!file.exists(fname)) stop("File not found: ", fname)
  df <- read.table(fname, header = TRUE, sep = "\t", quote = "",
                   stringsAsFactors = FALSE, check.names = TRUE)
  df$Sample <- sample
  df
}

all_data <- map(samples$sample_id, load_sample) |> bind_rows()
message(sprintf("Loaded %d / %d samples.", n_distinct(all_data$Sample), nrow(samples)))

# Guard: a gene symbol should appear once per sample in the -A table; if StringTie
# produced several rows (e.g. a novel locus sharing the name), sum them and say so.
gene_fpkm <- function(gene) {
  d <- all_data |>
    dplyr::filter(Gene.Name == gene) |>
    dplyr::group_by(Sample) |>
    dplyr::summarise(n_rows = dplyr::n(), FPKM = sum(FPKM), .groups = "drop")
  if (any(d$n_rows > 1))
    message(sprintf("NOTE: %s has >1 row in: %s (FPKM summed)", gene,
                    paste(d$Sample[d$n_rows > 1], collapse = ", ")))
  stopifnot(setequal(d$Sample, samples$sample_id))
  d |> dplyr::select(Sample, FPKM)
}

mycn  <- gene_fpkm("MYCN")  |> dplyr::rename(MYCN  = FPKM)
chek1 <- gene_fpkm("CHEK1") |> dplyr::rename(CHEK1 = FPKM)

# ── Table 1: MYCN_FPKM.txt (Fig 2B input), in samples.tsv order ──────────────
mycn_tbl <- samples |>
  dplyr::select(Sample = sample_id, MYCN_amp) |>
  dplyr::left_join(mycn, by = "Sample") |>
  dplyr::mutate(MYCN = round(MYCN, 2),
                MYCN_amp = MYCN_amp == "amp") |>
  dplyr::select(Sample, MYCN, MYCN_amp)
write_tsv(mycn_tbl, file.path(OUT_DIR, "MYCN_FPKM.txt"))

# ── Table 2: MYCN + CHEK1 FPKM per line (Fig 1A tracks) ──────────────────────
fpkm_tbl <- samples |>
  dplyr::select(Sample = sample_id, ecDNA_group) |>
  dplyr::left_join(mycn,  by = "Sample") |>
  dplyr::left_join(chek1, by = "Sample")
write_csv(fpkm_tbl, file.path(OUT_DIR, "MYCN_CHEK1_FPKM.csv"))

# ── Table 3: CHEK1_FPKM.csv (legacy name) ────────────────────────────────────
chek1_fpkm <- fpkm_tbl |>
  dplyr::select(Sample, ecDNA_group, FPKM = CHEK1) |>
  dplyr::arrange(dplyr::desc(ecDNA_group), Sample)   # ecDNA+ first, then ecDNA-
write_csv(chek1_fpkm, file.path(OUT_DIR, "CHEK1_FPKM.csv"))
print(chek1_fpkm, n = Inf)

cat(sprintf("\nMYCN FPKM  ecDNA+: median %.1f (range %.1f-%.1f) | ecDNA-: median %.1f (range %.1f-%.1f)\n",
            median(fpkm_tbl$MYCN[fpkm_tbl$ecDNA_group == "ecDNA+"]),
            min(fpkm_tbl$MYCN[fpkm_tbl$ecDNA_group == "ecDNA+"]),
            max(fpkm_tbl$MYCN[fpkm_tbl$ecDNA_group == "ecDNA+"]),
            median(fpkm_tbl$MYCN[fpkm_tbl$ecDNA_group == "ecDNA-"]),
            min(fpkm_tbl$MYCN[fpkm_tbl$ecDNA_group == "ecDNA-"]),
            max(fpkm_tbl$MYCN[fpkm_tbl$ecDNA_group == "ecDNA-"])))

# ── Optional: exploratory RS-gene heatmap (not a manuscript figure) ──────────
if (DRAW_HEATMAP) {
  rs_genes <- readLines(RS_GENES_TXT)
  rs_genes <- rs_genes[nzchar(rs_genes)]
  message(sprintf("RS panel: %d genes from %s", length(rs_genes), RS_GENES_TXT))

  rs_data <- all_data |>
    dplyr::filter(Gene.Name %in% rs_genes) |>
    dplyr::select(Sample, Gene.Name, FPKM) |>
    dplyr::mutate(log2FPKM = log2(FPKM + 1), ecDNA_group = ecdna_group[Sample])

  missing <- setdiff(rs_genes, unique(rs_data$Gene.Name))
  if (length(missing)) warning("Genes not found in data: ", paste(missing, collapse = ", "))
  stopifnot(!any(is.na(rs_data$ecDNA_group)))

  mat <- rs_data |>
    dplyr::select(Sample, Gene.Name, log2FPKM) |>
    tidyr::pivot_wider(names_from = Sample, values_from = log2FPKM,
                       values_fn = sum, values_fill = NA) |>
    tibble::column_to_rownames("Gene.Name") |>
    as.matrix()

  col_anno <- data.frame(ecDNA = ecdna_group[colnames(mat)], row.names = colnames(mat))
  anno_colors <- list(ecDNA = c("ecDNA+" = "#E64B35", "ecDNA-" = "#4DBBD5"))
  col_order <- order(col_anno$ecDNA, decreasing = TRUE)           # ecDNA+ first

  png(file.path(OUT_DIR, "RS_heatmap.png"), width = 10, height = 6, units = "in", res = 300)
  pheatmap(
    mat[, col_order],
    annotation_col    = col_anno[col_order, , drop = FALSE],
    annotation_colors = anno_colors,
    cluster_cols      = FALSE,
    cluster_rows      = TRUE,
    scale             = "row",
    color             = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
    border_color      = NA,
    fontsize          = 9,
    main              = sprintf("RS gene expression (n = %d; row z-score of log2[FPKM+1])", nrow(mat)),
    angle_col         = 45
  )
  dev.off()
  message("Saved: RS_heatmap.png (exploratory)")
}

message("Saved: MYCN_FPKM.txt, MYCN_CHEK1_FPKM.csv, CHEK1_FPKM.csv")
