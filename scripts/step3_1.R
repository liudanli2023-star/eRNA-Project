################################################################################
# Script: step3_1.R
# Purpose: Correct enhancer-gene correlation P values, calculate significance
#          scores, and filter high-confidence enhancer-gene candidate pairs.
# Inputs : 04.all_gene_enhancer_pair_pvalues.RData
#          04.all_gene_enhancer_pair_corrs.RData
#          02.GRCm39_enhancers_gene_1mb.txt
#          04.GRCm39_enhencer_noexon.bed
# Outputs: significant_pairs_all.RData
#          03.sig_enhancer_gene_raw_tpm.txt
#          filtered enhancer-gene pair tables under different thresholds
# Note   : The script uses chunk-wise FDR correction to reduce memory pressure.
################################################################################


gc()      # 清理内存

########################################step3#############################################
rm(list=ls())
getwd()
setwd("/home/dlliu/02_erna/13.enh_noncoding_file2")



# 设置并行线程数（根据CPU核心数调整，比如 8 核就设 plan(multisession, workers = 8)）
plan(multisession, workers = 64)


# 加载预计算的相关性和 p 值矩阵
load("01.file/04.all_gene_enhancer_pair_pvalues.RData")  # pvalueMatrix
load("01.file/04.all_gene_enhancer_pair_corrs.RData")    # corrMatrix

# 处理缺失值
pvalueMatrix[is.na(pvalueMatrix)] <- 1  # 将 NA p 值设为 1
corrMatrix[is.na(corrMatrix)] <- 0      # 将 NA 相关性设为 0

# FDR校正
#pvalueMatrix <- t(apply(pvalueMatrix, 1, p.adjust, method="fdr"))

## 方法3:FDR校正分步进行，减少内存，逐块处理 + 写入临时文件
library(future.apply)
plan(multisession, workers = 64)

n_rows <- nrow(pvalueMatrix)
chunk_size <- 1000
chunks <- split(1:n_rows, ceiling(seq_along(1:n_rows)/chunk_size))

# 创建结果文件
out_file <- "pvalueMatrix_corrected.bin"
file.create(out_file)

# 并行处理，每块写入临时文件
for (i in seq_along(chunks)) {
  cat("Processing chunk", i, "/", length(chunks), "\n")
  corrected <- t(apply(pvalueMatrix[chunks[[i]], , drop=FALSE], 1, p.adjust, method="fdr"))
  
  # 追加写入文件（用 readr::write_tsv 或 data.table::fwrite）
  save(corrected, file = paste0("chunk_", i, ".RData"))
}

# 合并所有块
corrected_chunks <- lapply(1:length(chunks), function(i) {
  load(paste0("chunk_", i, ".RData"))
  corrected
})
pvalueMatrix<- do.call(rbind, corrected_chunks)

## 可选 保存 pvalueMatrix
save(pvalueMatrix, file = "03.file/pvalueMatrix_FDR_corrected.RData")

##筛选显著增强子-基因对
# 筛选正相关（> 0.3）的增强子-基因对
row.col.indexes <- which(corrMatrix > 0.3, arr.ind = TRUE)

## 可选 保存 row.col.indexes**
save(row.col.indexes, file = "03.file/row_col_indexes_corr_gt0.3.RData")

# 获取增强子 ID 和基因 ID
enhancer_ids <- rownames(corrMatrix)[row.col.indexes[,1]]
gene_ids <- colnames(corrMatrix)[row.col.indexes[,2]]

# 提取相关性和 p 值
corr_values <- corrMatrix[row.col.indexes]
p_values <- pvalueMatrix[row.col.indexes]

# 计算显著性得分: 相关性 * (-log10(p值))
significance_scores <- corr_values * (-log10(p_values + 1e-100))

# 组合数据
significant_pairs <- data.frame(
  Enhancer = enhancer_ids,
  Gene = gene_ids,
  Correlation = corr_values,
  P_value = p_values,
  Significance_Score = significance_scores
)

## 可选 保存 significant_pairs（原始版）##
save(significant_pairs, file = "03.file/significant_pairs_raw.RData")
## write.table(significant_pairs, "03.file/significant_pairs_raw.txt",
            #row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")




######################################## restart from raw, 内存不够 #############################################
rm(list = ls())   # 清空内存
gc()              # 回收内存

setwd("/home/dlliu/02_erna/13.enh_noncoding_file2")

## 1. 加载已经保存好的 significant_pairs_raw
load("03.file/significant_pairs_raw.RData")   # 得到对象 significant_pairs


# 按得分排序，取top 5%
#significant_pairs <- significant_pairs[order(-significant_pairs$Significance_Score), ]
#top_n <- round(nrow(significant_pairs) * 0.05)
#significant_pairs <- significant_pairs[1:top_n, ]


# 去掉基因版本号
significant_pairs$Gene <- gsub("\\.\\d+", "", significant_pairs$Gene)

# 读取你的03.enh_gene.txt，筛选只保留这批
enh_gene <- read.table("02.file/02.GRCm39_enhancers_gene_1mb.txt", header=FALSE, stringsAsFactors=FALSE, col.names=c("Enhancer", "Gene"))
enh_gene$Gene <- gsub("\\.\\d+", "", enh_gene$Gene)

significant_pairs <- merge(significant_pairs, enh_gene, by=c("Enhancer", "Gene"))

## 可选 保存 significant_pairs（所有版）##
save(significant_pairs, file = "03.file/significant_pairs_all.RData")

# 合并增强子位置信息
enhancer_positions <- read.table("02.file/04.GRCm39_enhencer_noexon.bed", header=FALSE, stringsAsFactors=FALSE, 
                                 col.names=c("Enhancer", "Chr", "Start", "End"))
significant_pairs <- merge(significant_pairs, enhancer_positions, by="Enhancer", all.x=TRUE)

# 保存结果
write.table(significant_pairs, "03.sig_enhancer_gene_raw_tpm.txt",
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")



######################## 筛选条件：Correlation > 0.85 且 P_value < 0.01 ##########################
df <- read.table("03.sig_enhancer_gene_raw_tpm.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 筛选条件：Correlation > 0.85 且 P_value < 0.01
df_filtered <- subset(df, Correlation > 0.85 & P_value < 0.01)

# 保存筛选后的结果
write.table(df_filtered, "03.sig_enhancer_gene_tpm_filtered_0.85_0.01.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")



# 筛选条件：Correlation > 0.7 且 P_value < 0.05
df_filtered <- subset(df, Correlation > 0.7 & P_value < 0.05)

# 保存筛选后的结果
write.table(df_filtered, "03.sig_enhancer_gene_tpm_filtered_0.7_0.05.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")



# 筛选条件：Correlation > 0.7 且 P_value < 0.01
df_filtered <- subset(df, Correlation > 0.7 & P_value < 0.01)

# 保存筛选后的结果
write.table(df_filtered, "03.sig_enhancer_gene_tpm_filtered_0.7_0.01.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

