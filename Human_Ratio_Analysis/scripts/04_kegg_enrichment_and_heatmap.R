# 04_kegg_enrichment_and_heatmap.R
#
# Purpose: KEGG enrichment and summary heatmap
# Dataset: Human brain nuclear/cytoplasmic proteomics ratio data
# Notes:
#   - Run from the repository root directory.
#   - Input data and intermediate cluster gene lists are not included.
#   - See README.md for expected input files and workflow.
#
# 安装 KEGG 依赖

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# === Step 1: 读取所有 cluster 的基因并转换为 ENTREZ ID ===
kegg_list <- list()

for (i in 1:12) {
  genes <- readLines(paste0("Cluster_", i, "_genes.txt"))
  genes <- unique(genes[genes != ""])
  
  # 转换 symbol 为 EntrezID
  gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_ids <- gene_df$ENTREZID
  
  # KEGG 富集
  ekegg <- enrichKEGG(gene = entrez_ids, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  if (!is.null(ekegg) && nrow(ekegg) > 0) {
    kegg_list[[paste0("Cluster_", i)]] <- ekegg
  }
}

# === Step 2: 可视化 top10 term（每个 cluster 单独） ===
for (name in names(kegg_list)) {
  ekegg <- kegg_list[[name]]
  p <- barplot(ekegg, showCategory = 10, title = name, font.size = 12) +
    theme_minimal(base_size = 14)
  ggsave(paste0(name, "_KEGG_Top10.pdf"), p, width = 8, height = 6)
}


# 合并所有 KEGG 结果
kegg_all <- dplyr::bind_rows(
  lapply(names(kegg_list), function(cl) {
    df <- as.data.frame(kegg_list[[cl]])
    df$Cluster <- cl
    return(df)
  })
)

# 保存
write.csv(kegg_all, "KEGG_Merged_Clusters.csv", row.names = FALSE)

################################################绘制汇总KEGG热图

packages <- c("tidyverse", "ggplot2", "RColorBrewer", "scales", "reshape2")
lapply(packages, require, character.only = TRUE)

# === Step 2: 读取合并后的 KEGG 富集表 ===
kegg <- read.csv("KEGG_Merged_Clusters.csv", stringsAsFactors = FALSE)

# === Step 3: 为每个 cluster 选择 top5 term ===
kegg_top5 <- kegg %>%
  group_by(Cluster) %>%
  slice_min(order_by = p.adjust, n = 5) %>%
  ungroup()

# === Step 4: 构建热图矩阵（-log10FDR）===
kegg_top5 <- kegg_top5 %>%
  mutate(log10FDR = -log10(p.adjust + 1e-10))

# 设置 cluster 顺序
kegg_top5$Cluster <- factor(kegg_top5$Cluster,
                            levels = paste0("Cluster_", 1:10))

# 构建热图矩阵（行为 pathway，列为 cluster）
mat <- kegg_top5 %>%
  dplyr::select(Description, Cluster, log10FDR) %>%
  pivot_wider(names_from = Cluster, values_from = log10FDR, values_fill = 0)

mat_df <- as.data.frame(mat)
rownames(mat_df) <- mat_df$Description
mat_df <- mat_df[ , -1]  # 去掉 Description 列
desired_order <- paste0("Cluster_", 1:12)
mat_df <- mat_df[, desired_order[desired_order %in% colnames(mat_df)]]

# === Step 6: 标记关键 pathway（红色）===
highlight_terms <- c("Nucleocytoplasmic transport", "Spliceosome", "mRNA surveillance pathway")
highlight_rows <- rownames(mat_df) %>%
  str_to_lower() %>%
  str_detect(paste(tolower(highlight_terms), collapse = "|"))

# === Step 7: 绘图（热图，1:1，低饱和色）===
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
col_fun <- colorRamp2(c(0, max(mat_df)), c("#f2f2f2", "#4575b4"))

ht <- Heatmap(
  mat_df,
  name = "-log10(FDR)",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8, col = ifelse(highlight_rows, "red", "black")),
  column_names_gp = gpar(fontsize = 11),
  heatmap_width = unit(10, "cm"),
  heatmap_height = unit(10, "cm"),
  row_title = "KEGG Pathway",
  column_title = "Cluster",
  rect_gp = gpar(col = "grey80", lwd = 0.5)
)

pdf("KEGG_Top5_Heatmap_R.pdf", width = 7, height = 7)
draw(ht)
dev.off()

#######################20251106-KEGG聚类修改
# === Step 2: 读取合并后的 KEGG 富集表 ===
# === Step 1: 读取 KEGG 富集文件 ===

library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)

kegg <- read.csv("KEGG_Merged_Clusters.csv", stringsAsFactors = FALSE)

# === Step 2: 修复缺失或空行 ===

# 填充空 Description

kegg <- kegg %>%
  mutate(
    Description = ifelse(is.na(Description) | Description == "", NA, Description),
    p.adjust = ifelse(is.na(p.adjust) | p.adjust == "", 1, p.adjust)
  )

# === Step 3: 取每个 cluster 的 top5 term ===

kegg_top5 <- kegg %>%
  group_by(Cluster) %>%
  slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
  ungroup()

# === Step 4: 计算 -log10(FDR) ===

kegg_top5 <- kegg_top5 %>%
  mutate(log10FDR = -log10(as.numeric(p.adjust) + 1e-10))

# === Step 5: pivot_wider，自动处理重复项 ===

mat <- kegg_top5 %>%
  dplyr::select(Description, Cluster, log10FDR) %>%
  pivot_wider(
    names_from = Cluster,
    values_from = log10FDR,
    values_fill = list(log10FDR = 0),
    values_fn = list(log10FDR = max)
  )

# === Step 6: 生成矩阵并补齐缺失 cluster ===

mat_df <- as.data.frame(mat)
mat_df <- na.omit(mat_df)
rownames(mat_df) <- mat_df$Description
mat_df <- mat_df[, -1]

desired_order <- paste0("Cluster_", 1:12)
missing_clusters <- setdiff(desired_order, colnames(mat_df))

# 补齐缺失的 cluster 列

if (length(missing_clusters) > 0) {
  for (mc in missing_clusters) mat_df[[mc]] <- 0
}
mat_df <- mat_df[, desired_order]

# === Step 7: 标记关键通路 ===

highlight_terms <- c("Nucleocytoplasmic transport", "Spliceosome", "mRNA surveillance pathway")
highlight_rows <- rownames(mat_df) %>%
  str_to_lower() %>%
  str_detect(paste(tolower(highlight_terms), collapse = "|"))

# === Step 8: 绘图 ===

col_fun <- colorRamp2(c(0, max(mat_df)), c("#f2f2f2", "#4575b4"))

ht <- Heatmap(
  mat_df,
  name = "-log10(FDR)",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8, col = ifelse(highlight_rows, "red", "black")),
  column_names_gp = gpar(fontsize = 11),
  heatmap_width = unit(10, "cm"),
  heatmap_height = unit(10, "cm"),
  row_title = "KEGG Pathway",
  column_title = "Cluster",
  rect_gp = gpar(col = "grey80", lwd = 0.5)
)

pdf("KEGG_Top5_Heatmap_R_fixed.pdf", width = 7, height = 7)
draw(ht)
dev.off()
