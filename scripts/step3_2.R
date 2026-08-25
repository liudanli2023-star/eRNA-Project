################################################################################
# Script: step3_2.R
# Purpose: Annotate each significant enhancer-gene pair with the developmental
#          stage at which the enhancer shows its highest TPM value. The script
#          also summarizes one-enhancer/multiple-gene relationships.
# Inputs : 04.Mus_Brain_e10.5_9w_enhancer_TPM_0.1.txt
#          03.sig_enhancer_gene_tpm_filtered_0.85_0.01.txt
# Outputs: 04.sig_enhancer_gene_stage.txt
#          04.sig_enhancer_mergeGene_stage.txt
################################################################################

rm(list=ls())
getwd()
setwd("/home/dlliu/02_erna/13.enh_noncoding_file2/")


# 1. 加载必要库
tools <- c("dplyr")
new.packages <- tools[!(tools %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(dplyr)

# 1. 读取增强子表达矩阵
enhancer_expr <- read.table(
  "01.file/04.Mus_Brain_e10.5_9w_enhancer_TPM_0.1.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 2. 读取显著增强子-基因对
sig_pairs <- read.table(
  "03.file/03.sig_enhancer_gene_tpm_filtered_0.85_0.01.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 3. 提取显著增强子表达矩阵
expr_subset <- enhancer_expr[rownames(enhancer_expr) %in% sig_pairs$Enhancer, ]
expr_mat <- apply(expr_subset, 2, as.numeric)
rownames(expr_mat) <- rownames(expr_subset)

# 4. 计算每个增强子表达最高的阶段
max_idx <- apply(expr_mat, 1, which.max)
stages <- colnames(expr_mat)[max_idx]

enhancer.stages <- data.frame(
  Enhancer = rownames(expr_mat),
  Stage = stages,
  stringsAsFactors = FALSE
)

# 5. 读取原始增强子-基因关系
netData <- read.table(
  file = "03.file/03.sig_enhancer_gene_tpm_filtered_0.85_0.01.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 6. 添加阶段信息
matchIdx <- match(netData$Enhancer, enhancer.stages$Enhancer)
netData$Stage <- enhancer.stages$Stage[matchIdx]

# 7. 合并后直接保存总文件
write.table(
  netData,
  file = "03.file/04.sig_enhancer_gene_stage.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)




############################################
###### 数据整理,同一enhancer对应多个gene
library(dplyr)
library(readr)

# 读取数据
data <- read.table("03.file/04.sig_enhancer_gene_stage.txt", header = TRUE, sep = "\t")

# 合并基因：按 Enhancer + Stage 分组，基因合并，其他列取 unique（因为同一 enhancer 多基因时它们相同）
merged_data <- data %>%
  group_by(Enhancer, Stage, Chr, Start, End) %>%
  summarise(
    Gene = paste(unique(Gene), collapse = ","),
    Correlation = paste(Correlation, collapse = ","),
    P_value = paste(P_value, collapse = ","),
    Significance_Score = paste(Significance_Score, collapse = ","),
    .groups = "drop"
  )

# 保存结果
write.table(merged_data, "04.sig_enhancer_mergeGene_stage.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# write_csv(merged_data, "03.sig_enhancer_mergeGene_stage.csv")
