# Usage: Rscript 05_run_shatterseek.R <sample> <sv_dir> <cn_dir> <results_dir> [chromosomes]
#   <sv_dir>      folder with sv_<sample>.txt   (03_format_sv_for_shatterseek.py)
#   <cn_dir>      folder with <sample>_cn.txt   (04_format_cn_for_shatterseek.py)
#   <results_dir> receives <sample>_summary.txt (ShatterSeek chromSummary, read by 06_call_high_confidence_chromothripsis.ipynb)
#                 and <sample>_CT.pdf (per-chromosome plots)
#   [chromosomes] comma-separated chromosomes to plot; default "1,2,7,3,12" (the set used for the published run)
# Run once per sample (env ASCAT: R 4.2.2, ShatterSeek 1.1).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: Rscript 05_run_shatterseek.R <sample> <sv_dir> <cn_dir> <results_dir> [chromosomes]")
SAMPLE      <- args[1]
SV_DIR      <- file.path(args[2], "sv_")
CN_DIR      <- paste0(sub("/$", "", args[3]), "/")
Results_DIR <- paste0(sub("/$", "", args[4]), "/")
PLOT_CHRS   <- if (length(args) >= 5) trimws(strsplit(args[5], ",")[[1]]) else c("1", "2", "7", "3", "12")
dir.create(Results_DIR, showWarnings = FALSE, recursive = TRUE)
for (f in c(paste0(SV_DIR, SAMPLE, ".txt"), paste0(CN_DIR, SAMPLE, "_cn.txt")))
  if (!file.exists(f)) stop("missing input: ", f)

suppressPackageStartupMessages({
  library(ShatterSeek)
  library(gridExtra)
  library(cowplot)
  library(ggplot2)
  library(GenomicRanges)
})

sv = read.delim(paste0(SV_DIR, SAMPLE,'.txt'), sep='\t')
sv['strand1'] = ifelse(sv$ct == '3to3', '+', ifelse(sv$ct=='5to5', '-', ifelse(sv$ct=='3to5', '+','-')))
sv['strand2'] = ifelse(sv$ct == '3to3', '+', ifelse(sv$ct=='5to5', '-', ifelse(sv$ct=='3to5', '-','+')))

SV_data <- SVs(chrom1=as.character(sv$chr),
                           pos1=as.numeric(sv$start),
                           chrom2=as.character(sv$chr2),
                           pos2=as.numeric(sv$end),
                           SVtype=as.character(sv$svtype),
                           strand1=as.character(sv$strand1),
                           strand2=as.character(sv$strand2))

cna = read.delim(paste0(CN_DIR, SAMPLE,'_cn.txt'), sep='\t')

# Merge two adjacent regions that have the same copy number value
dd <- cna
dd$CN[dd$CN == 0] <- 150000
dd$CN[is.na(dd$CN)] <- 0
dd <- as(dd,"GRanges")
cov <- coverage(dd,weight = dd$CN)
dd1 <- as(cov,"GRanges")
dd1 <- as.data.frame(dd1)
dd1 <- dd1[dd1$score !=0,]
dd1 = dd1[,c(1,2,3,6)]
names(dd1) <- names(cna)[1:4]
dd1$CN[dd1$CN == 150000] <- 0
cna= dd1; rm(dd)

CN_data <- CNVsegs(chrom=as.character(cna$chromosome),
                                   start=cna$start,
                                   end=cna$end,
                                   total_cn=cna$CN)


chromothripsis <- shatterseek(
                SV.sample=SV_data,
                seg.sample=CN_data,
                genome="hg38")

file_path <- paste0(Results_DIR, SAMPLE, "_summary.txt")
write.table(chromothripsis@chromSummary, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)

# One 4-row panel (SV arcs, CN, oscillation states, breakpoints) per chromosome in PLOT_CHRS
pdf(NULL)                                   # arrangeGrob would otherwise open Rplots.pdf in the working directory
plot_list <- lapply(PLOT_CHRS, function(chr) {
  p <- plot_chromothripsis(ShatterSeek_output = chromothripsis, chr = chr, sample_name = SAMPLE, genome = "hg38")
  arrangeGrob(p[[1]], p[[2]], p[[3]], p[[4]], nrow = 4, ncol = 1, heights = c(0.2, .4, .4, .4))
})

pdf_path <- paste0(Results_DIR, SAMPLE, '_CT.pdf')

CT <- plot_grid(plotlist = plot_list)
ggsave(pdf_path, CT, device = "pdf", width = 10, height = 10)
invisible(dev.off())

