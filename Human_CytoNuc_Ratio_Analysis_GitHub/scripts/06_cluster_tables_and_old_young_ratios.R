# 06_cluster_tables_and_old_young_ratios.R
#
# Purpose: Combined cluster tables and selected old/young expression ratios
# Dataset: Human brain nuclear/cytoplasmic proteomics ratio data
# Notes:
#   - Run from the repository root directory.
#   - Input data and intermediate cluster gene lists are not included.
#   - See README.md for expected input files and workflow.
#
# 初始化一个空列表来存储所有数据框
all_data <- list()
# 循环处理每个文件
for (i in 1:12) {
  filename <- paste0("Cluster_", i, "_genes.txt")
  
  # 检查文件是否存在
  if (file.exists(filename)) {
    # 读取文件，假设文件只有一列基因名
    cluster_data <- read.table(filename, header = FALSE, stringsAsFactors = FALSE)
    colnames(cluster_data) <- "gene"
    
    # 添加cluster列
    cluster_data$cluster <- paste0("Cluster_", i)
    
    # 添加到列表
    all_data[[i]] <- cluster_data
  } else {
    warning(paste("File not found:", filename))
  }
}
# 合并所有数据框
combined_df <- do.call(rbind, all_data)
# 查看结果
head(combined_df)

write.csv(combined_df,file = "combined_cluster_Human.csv",row.names = F)


#########################################################
#按照cluster合并原始数据
# 加载需要的包
library(readr)
library(dplyr)
library(tibble) 
# 读取原始数据（按基因为行，13个年龄组为列）
original_data <- read_csv("Human_data_grouped_13ages.csv")

# 设置行名为gene symbol（假设第一列是gene名）
class(original_data)
original_data <- original_data %>%
  column_to_rownames(var = colnames(original_data)[1])

# 定义读取基因列表的函数
read_gene_list <- function(filepath) {
  scan(filepath, what = character(), sep = "\n", quiet = TRUE)
}

# 分别读取Cluster的基因列表
cluster6_genes <- read_gene_list("Cluster_6_genes.txt")
cluster9_genes <- read_gene_list("Cluster_9_genes.txt")
cluster10_genes <- read_gene_list("Cluster_10_genes.txt")

# 取交集后提取对应数据（注意检查基因名是否完全一致）
cluster6_data <- original_data[intersect(cluster6_genes, rownames(original_data)), ]
cluster9_data <- original_data[intersect(cluster9_genes, rownames(original_data)), ]
cluster10_data <- original_data[intersect(cluster10_genes, rownames(original_data)), ]


# 设置年龄组
#Young (20~29)
#Old (70~79)
young_ages <- c("21", "23", "25")
old_ages <- c("70", "73")

# 定义函数：计算old/young比值并提取top10,###排除0值
# 定义函数：仅在表达 > 0 的样本上计算平均
get_top10_ratio_genes <- function(cluster_data) {
  cluster_data <- cluster_data %>%
    mutate(
      Young_group_mean = rowMeans(select(., all_of(young_ages)), na.rm = TRUE),
      Old_group_mean = rowMeans(select(., all_of(old_ages)), na.rm = TRUE)
    ) %>%
    filter(Young_group_mean > 0) %>%  # 只保留年轻组有表达的行
    mutate(Ratio_old_vs_young = Old_group_mean / Young_group_mean) %>%
    arrange(desc(Ratio_old_vs_young))
  
  return(cluster_data)
}


cluster6_data_extended <- get_top10_ratio_genes(cluster6_data)
cluster9_data_extended <- get_top10_ratio_genes(cluster9_data)
cluster10_data_extended <- get_top10_ratio_genes(cluster10_data)



write.csv(cluster6_data_extended, "cluster6_data_ratio.csv")
write.csv(cluster9_data_extended, "cluster9_data_ratio.csv")
write.csv(cluster10_data_extended, "cluster10_data_ratio.csv")
