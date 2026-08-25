################################################################################
# Script: step3_3.R
# Purpose: Add stage-matched enhancer TPM and gene TPM values to significant
#          enhancer-gene pairs.
# Inputs : 04.sig_enhancer_gene_stage.txt
#          04.Mus_Brain_e10.5_9w_enhancer_TPM_0.1.txt
#          04.mus_brain_counts_tpm_filter4.txt
# Outputs: 04.sig_enhancerTPM_geneTPM_stage.txt
#          03.sig_enhancerTPM_geneTPM_stage.csv
################################################################################

rm(list=ls())
getwd()
setwd("/home/dlliu/02_erna/13.enh_noncoding_file2/")

# 读取输入文件
sig_pairs <- read.table("03.file/04.sig_enhancer_gene_stage.txt", header = TRUE, sep = "\t", check.names = FALSE)
enhancer_expr <- read.table("01.file/04.Mus_Brain_e10.5_9w_enhancer_TPM_0.1.txt", header = TRUE, row.names = 1, check.names = FALSE)
gene_expr <- read.table("01.file/04.mus_brain_counts_tpm_filter4.txt", header = TRUE, row.names = 1, check.names = FALSE)

# 提取14个时间点
timepoints <- c("X0d", "e10.5", "e11.5", "e12.5", "e13.5", "e14.5",
                "e15.5", "e16.5", "e17.5", "e18.5", "X2w", "X3d", "X4w", "X9w")
enhancer_expr <- enhancer_expr[, timepoints]
gene_expr <- gene_expr[, timepoints]

# 定义函数：根据 Stage 提取对应 TPM
get_tpm <- function(enhancer_id, gene_id, stage) {
  enh_tpm <- if (enhancer_id %in% rownames(enhancer_expr)) enhancer_expr[enhancer_id, stage] else NA
  gene_tpm <- if (gene_id %in% rownames(gene_expr)) gene_expr[gene_id, stage] else NA
  return(c(enh_tpm, gene_tpm))
}

# 对每一行 sig_pairs 添加 enhancer_TPM 和 gene_TPM
tpm_values <- mapply(get_tpm, 
                     sig_pairs$Enhancer, 
                     sig_pairs$Gene, 
                     sig_pairs$Stage)

# 转置，合并结果
tpm_values <- t(tpm_values)
colnames(tpm_values) <- c("Enhancer_TPM", "Gene_TPM")

sig_pairs_out <- cbind(sig_pairs, tpm_values)

# 保存结果
write.table(sig_pairs_out,
            file = "03.file/04.sig_enhancerTPM_geneTPM_stage.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(sig_pairs_out,
            file = "03.file/03.sig_enhancerTPM_geneTPM_stage.csv",
            sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
