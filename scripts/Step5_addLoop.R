################################################################################
# Script: Step5_addLoop.R
# Purpose: Add Hi-C loop evidence to marker gene-associated enhancer-gene pairs.
#          Full_Match indicates the same enhancer-gene pair is supported by a
#          loop. Enhancer_Same_Gene_Different indicates the enhancer is present
#          in loop data but linked to a different gene.
# Inputs : 03.sig_enhancerTPM_geneTPM_stage_markerGene_filter.txt
#          03.mus_8w_loops_PE_pairs.txt
#          03.mus_e12.5_loops_PE_pairs.txt
#          Mus_musculus.GRCm39.113.chr.gtf
# Outputs: 04.sig_enhancerTPM_geneTPM_stage_markerGene3_Loop.txt/csv
#          04.sig_enhancerTPM_geneTPM_stage_markerGene_GenePos.txt/csv
################################################################################


rm(list=ls())
getwd()
setwd("/home/dlliu/02_erna/13.enh_noncoding_file2/")


# 指定多个 loops 文件路径
loop_files <- c("4.hic_file/03.mus_8w_loops_PE_pairs.txt", 
                "4.hic_file/03.mus_e12.5_loops_PE_pairs.txt")

# 读取并合并所有 loops 文件
all_loops <- do.call(rbind, lapply(loop_files, function(f) {
  df <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("LoopID", "Chr1", "Start1", "End1", "GeneID", "Gene", 
                    "Chr2", "Start2", "End2", "Enhancer")
  df$Source_File <- f   # 记录来源文件
  return(df)
}))

# 读取 sig 文件
sig <- read.table("03.file/03.sig_enhancerTPM_geneTPM_stage_markerGene_filter.txt", 
                  header = TRUE, stringsAsFactors = FALSE)

readLines("03.file/03.sig_enhancerTPM_geneTPM_stage_markerGene_filter.txt")[188]
strsplit(readLines("03.file/03.sig_enhancerTPM_geneTPM_stage_markerGene_filter.txt")[188], "\t")[[1]] |> length()

# 有上述报错，换一个读取方式
library(data.table)
sig <- fread("03.file/03.sig_enhancerTPM_geneTPM_stage_markerGene_filter.txt", fill = TRUE)


# -------------------------
# 初始化标注列
# -------------------------
sig$Match_Status <- "No_Match"
sig$Loop_Source <- NA   # 记录匹配到的 loops 文件来源

# -------------------------
# 1. 完全匹配：Enhancer + Gene
# -------------------------
full_match_idx <- paste(sig$Enhancer, sig$Gene) %in% paste(all_loops$Enhancer, all_loops$Gene)
sig$Match_Status[full_match_idx] <- "Full_Match"

# 标注来源文件（取第一个匹配到的来源）
for (i in which(full_match_idx)) {
  enh <- sig$Enhancer[i]
  gene <- sig$Gene[i]
  matched_files <- unique(all_loops$Source_File[all_loops$Enhancer == enh & all_loops$Gene == gene])
  sig$Loop_Source[i] <- paste(matched_files, collapse = ",")
}

# -------------------------
# 2. Enhancer 相同但 Gene 不同
# -------------------------
common_enhancers <- intersect(sig$Enhancer, all_loops$Enhancer)

for (enh in common_enhancers) {
  sig_genes <- unique(sig$Gene[sig$Enhancer == enh])
  loop_genes <- unique(all_loops$Gene[all_loops$Enhancer == enh])
  
  # 如果 enhancer 存在但 gene 不一致（不重叠），标记
  diff_genes <- setdiff(sig_genes, loop_genes)
  if (length(diff_genes) > 0) {
    sig$Match_Status[sig$Enhancer == enh & sig$Gene %in% diff_genes & sig$Match_Status == "No_Match"] <- "Enhancer_Same_Gene_Different"
    sig$Loop_Source[sig$Enhancer == enh & sig$Gene %in% diff_genes] <- 
      paste(unique(all_loops$Source_File[all_loops$Enhancer == enh]), collapse = ",")
  }
}

# -------------------------
# 输出结果
# -------------------------
write.table(sig, "04.sig_enhancerTPM_geneTPM_stage_markerGene3_Loop.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sig, "04.sig_enhancerTPM_geneTPM_stage_markerGene3_Loop.csv", 
            sep = ",", quote = FALSE, row.names = FALSE)



############################ 添加gene的chr chr_start, chr_end #######################
library(rtracklayer)
library(dplyr)

# 1. 读入 sig 文件
sig <- read.table("04.sig_enhancerTPM_geneTPM_stage_markerGene3_Loop.txt",
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2. 读入 GTF 注释文件
gtf <- rtracklayer::import("Mus_musculus.GRCm39.113.chr.gtf")
genes <- gtf[gtf$type == "gene"]

# 转换成 data.frame
genes_df <- data.frame(
  Gene = ifelse(!is.na(mcols(genes)$gene_name), mcols(genes)$gene_name, mcols(genes)$gene_id),
  Gene_Chr = as.character(seqnames(genes)),
  Gene_Start = start(genes),
  Gene_End = end(genes),
  stringsAsFactors = FALSE
)

# 3. 合并到 sig
sig_annot <- sig %>%
  left_join(genes_df, by = "Gene")

# 4. 调整列顺序：在 Gene 后插入 Gene_Chr, Gene_Start, Gene_End
col_order <- c("Enhancer", "Gene", "Gene_Chr", "Gene_Start", "Gene_End", colnames(sig)[3:ncol(sig)])
sig_final <- sig_annot[, col_order]

# 5. 输出结果
write.table(sig_final, "04.sig_enhancerTPM_geneTPM_stage_markerGene_GenePos.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sig_final, "04.sig_enhancerTPM_geneTPM_stage_markerGene_GenePos.csv",
            sep = ",", quote = FALSE, row.names = FALSE)
