# 01_human_ratio_trajectory_clustering.R
#
# Purpose: Human ratio preprocessing, hierarchical clustering, cluster-wise trajectory visualization
# Dataset: Human brain nuclear/cytoplasmic proteomics ratio data
# Notes:
#   - Run from the repository root directory.
#   - Input data and intermediate cluster gene lists are not included.
#   - See README.md for expected input files and workflow.
#
# === Step 1: 加载依赖包 ===
packages <- c("tidyverse", "zoo", "dynamicTreeCut", "ggplot2")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])
invisible(lapply(packages, library, character.only = TRUE))

# === Step 2: 读取数据 ===
data_raw <- read.csv("Human_data_grouped_13ages.csv", row.names = 1, check.names = FALSE)

# 确保年龄列排序
data_raw <- data_raw[, order(as.numeric(colnames(data_raw)))]

# === Step 3: 缺失值线性插值填补 ===
data_filled <- t(apply(data_raw, 1, function(x) {
  if (sum(!is.na(x)) >= 2) {
    zoo::na.approx(x, rule = 2)
  } else {
    rep(NA, length(x))  # 保证维度一致
  }
}))
# 修复列名
colnames(data_filled) <- colnames(data_raw)
rownames(data_filled) <- rownames(data_raw)
head(data_filled)[,1:6]
sum(is.na(data_filled))  # 检查还有多少 NA

# 删除仍含 NA 的蛋白
data_filled <- data_filled[rowSums(is.na(data_filled)) == 0, ]

# 检查剩余蛋白数量
nrow(data_filled)
# 先移除标准差为 0 的蛋白
row_sd <- apply(data_filled, 1, sd)
data_filled <- data_filled[row_sd != 0, ]
sum(row_sd == 0)  # 有多少蛋白是恒定值

# === Step 4: Z-score 标准化（按行） ===
data_scaled <- t(scale(t(data_filled)))

# === Step 5: 聚类分析 ===
dist_matrix <- dist(data_scaled, method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")
###自动聚类是25个cluster
#clusters <- dynamicTreeCut::cutreeDynamic(hc, distM = as.matrix(dist_matrix), verbose = 0)
#table(clusters)#########dynamicTreeCut::cutreeDynamic自动聚类
###########层次聚类
library(dendextend)
# 构建 dendrogram
dend <- as.dendrogram(hc)

#####切换为12个
# 替换 dynamicTreeCut 为 12 个 cluster 的硬切分
clusters <- cutree(hc, k = 12)

# 替换 0 为 NA（为了避免出现在图中）
clusters_filtered <- clusters
clusters_filtered[clusters_filtered == 0] <- NA
clusters_filtered <- factor(clusters_filtered)
table(clusters_filtered)
# 创建颜色映射表
cluster_colors <- RColorBrewer::brewer.pal(n = length(unique(na.omit(clusters_filtered))), name = "Set3")

# 将颜色分配给分支
dend <- color_branches(dend, k = length(cluster_colors), col = cluster_colors)

# 绘图
pdf("Hierarchical Clustering Dendrogram of Proteins.pdf")
plot(dend, main = "Hierarchical Clustering Dendrogram of Proteins", ylab = "Height")
dev.off()

length(unique(df_long$Gene))
###################调整dend层面分析cluster顺序，使其服从于下游的轨迹分析顺序
######20260127
library(dendextend)
library(RColorBrewer)
# 替换 cluster==0 为 NA（未分类蛋白）
# 替换 dynamicTreeCut 为 13 个 cluster 的硬切分
clusters <- cutree(hc, k = 12)#####这是重点，必要

# 替换 0 为 NA（为了避免出现在图中）
clusters_filtered <- clusters
clusters_filtered[clusters_filtered == 0] <- NA
clusters_filtered <- factor(clusters_filtered)
cluster_gene_list <- split(names(clusters_filtered), clusters_filtered)
traj_order <- names(table(clusters_filtered))
traj_order
# 新编号（按轨迹分析顺序）
new_labels <- paste0("C", seq_along(traj_order))
names(new_labels) <- traj_order

new_labels
clusters_final <- new_labels[as.character(clusters)]
table(clusters_final)

cluster_gene_list_final <- split(
  names(clusters_final),
  clusters_final
)

dend <- as.dendrogram(hc)

# 不显示蛋白名
labels(dend) <- rep("", length(labels(dend)))

# 按 cluster 上色（用原始结构）
cluster_colors <- RColorBrewer::brewer.pal(n = length(unique(na.omit(clusters_filtered))), name = "Set3")


# 将颜色分配给分支
dend <- color_branches(dend, k = length(cluster_colors), col = cluster_colors)
dend_order <- order.dendrogram(dend)

# 每个 cluster 在树上的中心位置
cluster_x <- tapply(
  seq_along(dend_order),
  clusters[dend_order],
  mean
)

# 用“新编号”替换原来的数字
cluster_labels <- new_labels[names(cluster_x)]

cut_height <- max(hc$height) * 0.6

pdf("Dendrogram_aligned_to_Trajectory_Order_Mouse.pdf", width = 6, height = 6)

plot(dend,
     main = "Hierarchical Clustering (Aligned to Trajectory Order)",
     ylab = "Height",
     leaflab = "none")

text(
  x = cluster_x,
  y = rep(cut_height, length(cluster_x)),
  labels = cluster_labels,
  cex = 1.3,
  font = 2
)

dev.off()

###################################################################
# 转为 tidy 数据框用于 ggplot
df_long <- data_scaled %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  mutate(Cluster = clusters_filtered) %>%
  pivot_longer(-c(Gene, Cluster), names_to = "Age", values_to = "Zscore") %>%
  mutate(Age = as.numeric(Age)) %>%
  filter(!is.na(Cluster))  # 排除未分簇的蛋白


# 按 cluster 分组导出基因名
cluster_gene_list <- split(names(clusters_filtered), clusters_filtered)
for (k in names(cluster_gene_list)) {
  writeLines(cluster_gene_list[[k]], paste0("Cluster_", k, "_genes.txt"))
}

#######################
cluster_sizes <- df_long %>%
  group_by(Cluster) %>%
  summarise(N = n_distinct(Gene)) %>%
  ungroup()

# 合并回主数据框用于绘图注释
df_long <- df_long %>%
  left_join(cluster_sizes, by = "Cluster")

#########
# == 筛选 Cluster 1 的下行趋势基因 ==
# 提取 Cluster 1 的数据
cluster1_long <- df_long %>% filter(Cluster == 1)

# 为每个基因计算线性回归斜率
cluster1_slopes <- cluster1_long %>%
  group_by(Gene) %>%
  summarise(slope = coef(lm(Zscore ~ Age))[2]) %>%
  filter(slope < 0)

# 保留斜率为负的基因
df_long <- df_long %>%
  filter(!(Cluster == 1 & !(Gene %in% cluster1_slopes$Gene)))

# === 导出筛选后的 Cluster 1 基因列表 ===
writeLines(cluster1_slopes$Gene, "Cluster_1_genes_filtered.txt")

# === 更新 Cluster 1 的 N 为筛选后数量 ===
n_cluster1 <- cluster1_slopes %>% nrow()

df_long <- df_long %>%
  mutate(N = ifelse(Cluster == 1, n_cluster1, N))

head(df_long)

#######再运行一遍cluster_sizes
cluster_sizes <- df_long %>%
  group_by(Cluster) %>%
  summarise(N = n_distinct(Gene)) %>%
  ungroup()

# === Step 1: 计算平均轨迹 ===
df_mean <- df_long %>%
  group_by(Cluster, Age) %>%
  summarise(mean_Z = mean(Zscore), .groups = "drop")

# === Step 2: 提取线性回归斜率 ===
cluster_slopes <- df_mean %>%
  group_by(Cluster) %>%
  summarise(Slope = coef(lm(mean_Z ~ Age))[2]) %>%
  arrange(desc(Slope))

# 按照 slope 给 cluster 重新排序
cluster_ordered <- cluster_slopes %>%
  mutate(Cluster = factor(Cluster, levels = Cluster))

# 标记关键 cluster：升最快前1个，降最快前1个（你可调数量）
key_clusters <- cluster_slopes %>%
  top_n(1, Slope) %>%
  bind_rows(cluster_slopes %>% top_n(-1, Slope)) %>%
  pull(Cluster)

# 创建带标签的 Cluster 名
cluster_slopes <- cluster_slopes %>%
  mutate(ClusterLabel = ifelse(Cluster %in% key_clusters,
                               paste0("Cluster ", Cluster, " *"),
                               paste0("Cluster ", Cluster)))

# 合并标签回 df_mean 和 cluster_sizes
df_mean <- df_mean %>%
  left_join(cluster_slopes, by = "Cluster")

cluster_sizes <- cluster_sizes %>%
  left_join(cluster_slopes, by = "Cluster")

pdf("Protein Nuclear Ratio Trajectories by Cluster_1.pdf")
ggplot(df_mean, aes(x = Age, y = mean_Z)) +
  stat_smooth(method = "loess", se = TRUE, span = 0.75,
              color = "black", fill = "grey80", size = 4, level = 0.95) +
  facet_wrap(~ fct_inorder(ClusterLabel), ncol = 4) +
  scale_y_continuous(limits = c(-2, 2)) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 13),
    strip.background = element_rect(fill = "grey90", color = "black", size = 0.5)
  ) +
  labs(title = "Protein Nuclear Ratio Trajectories by Cluster (Sorted by Trend)",
       y = "Z-score", x = "Age") +
  geom_text(
    data = cluster_sizes,
    aes(x = Inf, y = Inf, label = paste0("n=", N)),
    hjust = 1.1, vjust = 1.5,
    inherit.aes = FALSE,
    size = 5
  )
dev.off()

write.csv(cluster_sizes,file = "human_cluster_size_slope.csv",row.names = F)

##################################################################
