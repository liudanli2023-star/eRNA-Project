################################################################################
# Script: step4_markergene.R
# Purpose: Identify enhancer-gene pairs associated with mouse brain developmental
#          marker genes, add cell-type labels, and filter by enhancer TPM and
#          enhancer length.
# Inputs : 04.sig_enhancerTPM_geneTPM_stage.txt
#          marker_celltype_gene.txt (optional general marker table)
# Outputs: 03.sig_enhancerTPM_geneTPM_stage_markerGene.txt
#          03.sig_enhancerTPM_geneTPM_stage_markerGene_filter.txt
#          03.sig_enhancerTPM_geneTPM_stage_marker.txt
################################################################################

rm(list=ls())
getwd()
setwd("/home/dlliu/02_erna/13.enh_noncoding_file2/")

# 加载数据
df <- read.table("03.file/04.sig_enhancerTPM_geneTPM_stage.txt",
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 定义 marker gene 和 cell-type 对应表
marker_list <- list(
  "Neuron-Excitatory"   = c("Slc17a7","Slc17a6","Slc17a8","Neurod2"),
  "Neuron-Inhibitory"   = c("Gad1","Gad2","Slc32a1","Dlx2"),
  "Neuron-MSN"          = c("Dpp6","Drd2","Adora2a"),
  "Neurons"             = c("Rbfox3","Tubb3"),
  "Neuroblasts"         = c("Ascl1","Dcx"),
  "Glioblasts"          = c("Sox10","Tnc","Egfr"),
  "Radialglia"          = c("Fabp7","Pax6"),
  "OPCs"                = c("Pdgfra","Olig2","Vcan"),
  "Oligodendrocyte"     = c("Mog","Mbp","Ptgds"),
  "Microglial_cell"     = c("Cx3cr1","C1qa","Aif1","Hexb"),
  "Astrocyte"           = c("Aldh1l1", "Slc1a3", "Gfap","Aqp4"),
  "Macrophage"          = c("Apoe","Mrc1"),
  "Endothelial"         = c("Cldn5"),
  "IPCS"                = c("Tbr2"),
  "Early_cortical_neurons" = c("Tbr1","Bcl11b"),
  "Late_cortical_neurons"  = c("Satb2","Cux1","Cux2"),
  "L2/3_IT"             = c("Calb1"),
  "L4/5_IT"             = c("Cux2","Rorb"),
  "L5PT_1"              = c("Pou3f1"),
  "L6_IT_2"             = c("Slc24a2"),
  "L6_CT_1"             = c("Syt6"),
  "L5IT"                = c("Ptn"),
  "Apical progenitors" = c("Sox2"),
  "Intermediate progenitors" = c("Eomes"),
  "Migrating neurons" = c("Nrp1"),
  "Corticothalamic PN" = c("Tle4"),
  "ependymocytes" = c("Foxj1")
)


# 展开为 data.frame 方便匹配
marker_df <- do.call(rbind, lapply(names(marker_list), function(ct) {
  data.frame(CellType = ct, Gene = marker_list[[ct]], stringsAsFactors = FALSE)
}))

# 筛选：只保留 Gene 在 marker gene 中的行
merged <- merge(df, marker_df, by = "Gene")

# 调整列顺序：原始列 + CellType 最后一列
merged <- merged[, c(colnames(df), "CellType")]

# 保存结果
out_file <- "03.file/03.sig_enhancerTPM_geneTPM_stage_markerGene.txt"
write.table(merged, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("筛选完成，结果已保存到:", out_file, "\n")





rm(list=ls())
###### 筛选 eRNA_TPM 大于1，长度 > 180, < 10000

# 输入输出文件
input_file <- "03.file/03.sig_enhancerTPM_geneTPM_stage_markerGene.txt"
output_txt <- "03.file/03.sig_enhancerTPM_geneTPM_stage_markerGene_filter.txt"


# 读入数据
df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 计算长度
df$Length <- df$End - df$Start

# 筛选条件：Enhancer_TPM > 1 且 长度 > 180, < 10000
filtered <- subset(df, Enhancer_TPM > 1 & Length > 180 & Length < 10000)

# 保存为制表符分隔的 txt
write.table(filtered, file = output_txt, sep = "\t", quote = FALSE, row.names = FALSE)


cat("✅ 筛选完成，结果已保存到:\n", output_txt, "\n")









# ==========================
# 通用脚本：挑选 marker gene 并加上 CellType
# ==========================

# 输入文件路径
input_file <- "03.file/03.sig_enhancerTPM_geneTPM_stage.txt"
marker_file <- "03.file/marker_celltype_gene.txt"
output_file <- "03.file/03.sig_enhancerTPM_geneTPM_stage_marker.txt"

# 读入主文件
df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 读入 marker 文件（要求两列：CellType, Gene）
marker_df <- read.table(marker_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 合并：只保留 Gene 在 marker gene 中的行
merged <- merge(df, marker_df, by = "Gene")

# 调整列顺序：原始列 + CellType 最后一列
merged <- merged[, c(colnames(df), "CellType")]

# 保存结果
write.table(merged, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ 筛选完成，结果已保存到:", output_file, "\n")
