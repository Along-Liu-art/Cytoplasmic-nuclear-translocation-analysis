# 03_cluster_amplitude_variability.R
#
# Purpose: Cluster amplitude and variability indices
# Dataset: Human brain nuclear/cytoplasmic proteomics ratio data
# Notes:
#   - Run from the repository root directory.
#   - Input data and intermediate cluster gene lists are not included.
#   - See README.md for expected input files and workflow.
#
# 加载数据
expr <- read.csv("Human_data_grouped_13ages.csv", row.names = 1, check.names = FALSE)
expr <- expr[complete.cases(expr), ]  # 去除有 NA 的蛋白

# 按年龄顺序排列列名
expr <- expr[, order(as.numeric(colnames(expr)))]
# 读取每个 cluster 的基因
cluster_list <- list()
for (i in 1:12) {
  genes <- readLines(paste0("Cluster_", i, "_genes.txt"))
  cluster_list[[as.character(i)]] <- intersect(genes, rownames(expr))
}
#Step 3：计算 Amplitude 和 Variability 指数
amplitude_index <- numeric()##代表该 cluster 中蛋白的表达变化强度
variability_index <- numeric()##代表该 cluster 中蛋白之间在每个时间点的表达一致性（越高表示差异大）

for (i in 1:12) {
  genes <- cluster_list[[as.character(i)]]
  mat <- expr[genes, , drop = FALSE]
  
  # Amplitude: 蛋白 max-min，再求平均
  amp <- apply(mat, 1, function(x) max(x) - min(x))
  amplitude_index[i] <- mean(amp)
  
  # Variability: 每个时间点上蛋白之间的 sd，再取平均
  var <- apply(mat, 2, sd)
  variability_index[i] <- mean(var)
}

df_index <- data.frame(
  Cluster = factor(1:12),
  Amplitude = amplitude_index,
  Variability = variability_index
)

####
library(ggplot2)
library(RColorBrewer)
pastel12 <- c(
  "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4",
  "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
  "#F2F2F2",
  "#D0E6FA",  # 新增灰蓝
  "#E6CCE6",  # 新增淡紫
  "#FFF2CC"   # 新增奶黄色
)

pdf("Cluster Stability vs Variability.pdf")
ggplot(df_index, aes(x = Amplitude, y = Variability, label = Cluster)) +
  geom_point(aes(color = Cluster), size = 6) +
  geom_text(color = "black", size = 4.5, fontface = "bold") +
  scale_color_manual(values = pastel12) +
  theme_classic(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 13)
  ) +
  labs(title = "Cluster Stability vs Variability",
       x = "Amplitude Index",
       y = "Variability Index") +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1)
dev.off()
