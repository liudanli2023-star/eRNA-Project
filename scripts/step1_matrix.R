################################################################################
# Script: step1_matrix.R
# Purpose: Calculate Pearson correlation and P-value matrices between enhancer
#          TPM profiles and gene TPM profiles across mouse brain developmental
#          stages.
# Inputs : 04.Mus_Brain_e10.5_9w_enhancer_TPM.txt
#          04.mus_brain_counts_tpm_filter4.txt
# Outputs: 04.all_gene_enhancer_pair_corrs.RData
#          04.all_gene_enhancer_pair_pvalues.RData
# Note   : Update setwd() and input filenames before running in a new environment.
################################################################################

rm(list = ls())

getwd()
setwd('/home/mantt/04.erna_file/')

library(data.table)
library(foreach)
library(doParallel)

# --------------------------
# 读取数据
# --------------------------
enhancer_df <- fread("04.Mus_Brain_e10.5_9w_enhancer_TPM.txt", header = TRUE)
gene_df     <- fread("04.mus_brain_counts_tpm_filter4.txt", header = TRUE)

# 转换为矩阵
enhancerMat <- as.matrix(enhancer_df[ , -1, with = FALSE])
rownames(enhancerMat) <- enhancer_df[[1]]

geneMat <- as.matrix(gene_df[ , -1, with = FALSE])
rownames(geneMat) <- gene_df[[1]]

cat("增强子矩阵维度:", dim(enhancerMat), "\n")
cat("基因矩阵维度:", dim(geneMat), "\n")

# --------------------------
# 分块参数
# --------------------------
block_size <- 2000
n_blocks   <- ceiling(nrow(enhancerMat) / block_size)
cat("总块数:", n_blocks, "\n")

n <- ncol(enhancerMat)
df <- n - 2

# --------------------------
# 并行设置
# --------------------------
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)


# --------------------------
# 并行分块计算
# --------------------------
results <- foreach(i = 1:n_blocks, .packages = "stats") %dopar% {
  start_idx <- (i - 1) * block_size + 1
  end_idx   <- min(i * block_size, nrow(enhancerMat))
  idx_range <- start_idx:end_idx
  
  enh_block <- enhancerMat[idx_range, , drop = FALSE]
  
  # 相关性
  corr_block <- cor(t(enh_block), t(geneMat), method = "pearson")
  
  # p 值
  t_values <- corr_block * sqrt(df / (1 - corr_block^2))
  p_block  <- 2 * pt(-abs(t_values), df)
  
  list(corr = corr_block, pval = p_block, enh_ids = rownames(enh_block))
}

stopCluster(cl)

# --------------------------
# 合并结果
# --------------------------
corrMatrix   <- do.call(rbind, lapply(results, function(x) x$corr))
pvalueMatrix <- do.call(rbind, lapply(results, function(x) x$pval))
rownames(corrMatrix)   <- unlist(lapply(results, function(x) x$enh_ids))
rownames(pvalueMatrix) <- unlist(lapply(results, function(x) x$enh_ids))
colnames(corrMatrix)   <- rownames(geneMat)
colnames(pvalueMatrix) <- rownames(geneMat)

# --------------------------
# 保存结果
# --------------------------
save(corrMatrix,   file = "04.all_gene_enhancer_pair_corrs.RData")
save(pvalueMatrix, file = "04.all_gene_enhancer_pair_pvalues.RData")

cat("✅ 并行分块完成！结果矩阵已保存到 RData 文件。\n")
