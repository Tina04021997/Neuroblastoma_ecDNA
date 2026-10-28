#!/usr/bin/env Rscript
# ============================================================
# Export StringTie FPKM tables used by the figures
# ============================================================
# Reads <sample>.gene.abundance.txt (StringTie -A output, 04_rnaseq/06_stringtie_batch{1,2}.sh) for the
# 21 lines listed in samples.tsv and writes:
#   MYCN_CHEK1_FPKM.csv  Sample / ecDNA_group / MYCN / CHEK1 -> 06_figures/00_build_figure_tables.py, which derives
#                                                             MYCN_FPKM.txt (Fig 2B), FPKM_CN_MYCN.txt (Fig 2A) and the
#                                                             MYCN_FPKM / CHEK1_FPKM columns of data/master.csv (Fig 1A)
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

# ── Paths: pass as arguments (or edit the defaults) ─────────────────────────
#   Rscript 03_export_FPKM_tables.R <stringtie_dir> <samples.tsv> <rs_genes.txt> <out_dir>
# <stringtie_dir> must hold the 21 <sample_id>.gene.abundance.txt files in ONE flat folder (06_stringtie_batch1.sh writes
# them flat; 06_stringtie_batch2.sh writes <sample>/<sample>.gene.abundance.txt -- collect both batches first).
args          <- commandArgs(trailingOnly = TRUE)
STRINGTIE_DIR <- if (length(args) >= 1) args[1] else "/path/to/stringtie_output"
SAMPLES_TSV   <- if (length(args) >= 2) args[2] else "/path/to/repo/samples.tsv"
RS_GENES_TXT  <- if (length(args) >= 3) args[3] else "/path/to/repo/data/deseq2_gsea/rs_genes.txt"   # written by 02_gsea_and_RS_genes.R
OUT_DIR       <- if (length(args) >= 4) args[4] else STRINGTIE_DIR
DRAW_HEATMAP  <- TRUE
if (any(grepl("^/path/to/", c(STRINGTIE_DIR, SAMPLES_TSV, RS_GENES_TXT))))
  stop("pass <stringtie_dir> <samples.tsv> <rs_genes.txt> <out_dir> as arguments (or edit the defaults in this script)")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
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

# (MYCN_FPKM.txt is no longer written here: 00_build_figure_tables.py derives it from MYCN_CHEK1_FPKM.csv, so the
#  repository has a single source for every figure table.)

# ── Table 1: MYCN + CHEK1 FPKM per line ──────────────────────────────────────
fpkm_tbl <- samples |>
  dplyr::select(Sample = sample_id, ecDNA_group) |>
  dplyr::left_join(mycn,  by = "Sample") |>
  dplyr::left_join(chek1, by = "Sample")
write_csv(fpkm_tbl, file.path(OUT_DIR, "MYCN_CHEK1_FPKM.csv"))

# ── Table 2: CHEK1_FPKM.csv (legacy name) ────────────────────────────────────
chek1_fpkm <- fpkm_tbl |>
  dplyr::select(Sample, ecDNA_group, FPKM = CHEK1) |>
  dplyr::arrange(dplyr::desc(ecDNA_group == "ecDNA+"), Sample)   # ecDNA+ first, then ecDNA- (locale-independent)
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
  stopifnot(!anyDuplicated(rs_genes), length(unique(rs_genes)) == 29L)   # audit point 9
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
  col_order <- order(col_anno$ecDNA == "ecDNA+", decreasing = TRUE)   # ecDNA+ first (locale-independent)

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

message("Saved: MYCN_CHEK1_FPKM.csv, CHEK1_FPKM.csv")
