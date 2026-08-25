################################################################################
# Script: step4_ASDgene.R
# Purpose: Screen conserved enhancer-gene pairs associated with autism spectrum
#          disorder (ASD)-related genes using the mouse SFARI gene list.
# Inputs : SFARI_gene_mouse.txt
#          01.sig_enhancerTPM_geneTPM_orthologs_raw.txt
# Output : 02.sig_enhancerTPM_geneTPM_ASDgene_mus.txt
# Note   : Gene symbols are case-sensitive and should be checked before matching.
################################################################################
rm(list=ls())
getwd()
setwd("/home/dlliu/02_erna/13.enh_noncoding_file2")



######################### 小鼠的ASD_gene
# 读入 SFARI gene list
sfari_mus_genes <- read.table("SFARI_gene_mouse.txt", header = FALSE, stringsAsFactors = FALSE)
sfari_mus_genes <- sfari_mus_genes$V1   # 提取基因名向量


# 读入 enhancer-gene 数据
orthologs <- read.table("01.sig_enhancerTPM_geneTPM_orthologs_raw.txt",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 筛选 Gene 列在 SFARI 基因列表中的行
mus_matched_rows <- orthologs[orthologs$Gene %in% sfari_mus_genes, ]


# 查看前几行
head(mus_matched_rows)

# 保存结果
write.table(mus_matched_rows, "02.sig_enhancerTPM_geneTPM_ASDgene_mus.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
