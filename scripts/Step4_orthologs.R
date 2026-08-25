################################################################################
# Script: Step4_orthologs.R
# Purpose: Screen mouse enhancer-gene pairs whose enhancer IDs overlap with
#          mouse-human orthologous enhancer coordinates.
# Inputs : 02.hg19_mus39_enhancer_overlap.bed
#          04.sig_enhancerTPM_geneTPM_stage.txt
# Output : 01.sig_enhancerTPM_geneTPM_orthologs_raw.txt
################################################################################


rm(list=ls())
getwd()
setwd("/home/dlliu/02_erna/13.enh_noncoding_file2/")


# 读入文件
bed_overlap <- read.table("5.orthologs_file/02.hg19_mus39_enhancer_overlap.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
sig_file <- read.table("03.file/04.sig_enhancerTPM_geneTPM_stage.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 取第8列（Enhancer ID）
overlap_ids <- bed_overlap[, 8]

# 筛选 03 文件中第一列匹配的行
matched_rows <- sig_file[sig_file$Enhancer %in% overlap_ids, ]

# 查看结果
head(matched_rows)

# 保存结果
write.table(matched_rows, "01.sig_enhancerTPM_geneTPM_orthologs_raw.txt", sep="\t", quote=FALSE, row.names=FALSE)
