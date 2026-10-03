library(ShatterSeek)
library('stringr')
library(gridExtra)
library(cowplot)
library('ggplot2')
library(GenomicRanges)

SAMPLE='NGP'
# ── Paths: edit before running ──────────────────────────────────────────────
SV_DIR='/path/to/shatterseek/sv_files/sv_'          # prefix of sv_<sample>.txt (03_format_sv_for_shatterseek.py)
CN_DIR='/path/to/shatterseek/cn_files/'             # <sample>_cn.txt (04_format_cn_for_shatterseek.py)
Results_DIR='/path/to/shatterseek/results/'

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

file_path <- paste0(SAMPLE, "_summary.txt")
# Write the data to a text file
write.table(chromothripsis@chromSummary, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)

### MODIFY THE FOLLOWING CHR2 TO YOUR TARGET ONES ###

plots_chr1 <- plot_chromothripsis(ShatterSeek_output = chromothripsis,
              chr = "1", sample_name=SAMPLE, genome="hg38")
plot_chr1 = arrangeGrob(plots_chr1[[1]],
                         plots_chr1[[2]],
                         plots_chr1[[3]],
                         plots_chr1[[4]],
                         nrow=4,ncol=1,heights=c(0.2,.4,.4,.4))

plots_chr2 <- plot_chromothripsis(ShatterSeek_output = chromothripsis,
              chr = "2", sample_name=SAMPLE, genome="hg38")
plot_chr2 = arrangeGrob(plots_chr2[[1]],
                         plots_chr2[[2]],
                         plots_chr2[[3]],
                         plots_chr2[[4]],
                         nrow=4,ncol=1,heights=c(0.2,.4,.4,.4))

plots_chr7 <- plot_chromothripsis(ShatterSeek_output = chromothripsis,
              chr = "7", sample_name=SAMPLE, genome="hg38")
plot_chr7 = arrangeGrob(plots_chr7[[1]],
                         plots_chr7[[2]],
                         plots_chr7[[3]],
                         plots_chr7[[4]],
                         nrow=4,ncol=1,heights=c(0.2,.4,.4,.4))

plots_chr3 <- plot_chromothripsis(ShatterSeek_output = chromothripsis,
              chr = "3", sample_name=SAMPLE, genome="hg38")
plot_chr3 = arrangeGrob(plots_chr3[[1]],
                         plots_chr3[[2]],
                         plots_chr3[[3]],
                         plots_chr3[[4]],
                         nrow=4,ncol=1,heights=c(0.2,.4,.4,.4))

plots_chr12 <- plot_chromothripsis(ShatterSeek_output = chromothripsis,
              chr = "12", sample_name=SAMPLE, genome="hg38")
plot_chr12 = arrangeGrob(plots_chr12[[1]],
                         plots_chr12[[2]],
                         plots_chr12[[3]],
                         plots_chr12[[4]],
                         nrow=4,ncol=1,heights=c(0.2,.4,.4,.4))


pdf_path <- paste0(Results_DIR, SAMPLE, '_CT.pdf')

CT <- plot_grid(plot_chr1, plot_chr2, plot_chr7, plot_chr3, plot_chr12)
ggsave(pdf_path, CT, device = "pdf", width = 10, height = 10)

