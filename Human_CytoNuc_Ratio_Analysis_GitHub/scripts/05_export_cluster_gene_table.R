# 05_export_cluster_gene_table.R
#
# Purpose: Export cluster membership into a long-format table
# Dataset: Human brain nuclear/cytoplasmic proteomics ratio data
# Notes:
#   - Run from the repository root directory.
#   - Input data and intermediate cluster gene lists are not included.
#   - See README.md for expected input files and workflow.
#
library(tidyverse)
# 加载数据
# 读取每个 cluster 的基因
cluster_list <- list()

for (i in 1:12) {
  
  genes <- readLines(
    paste0("Cluster_", i, "_genes.txt")
  )
  
  cluster_list[[as.character(i)]] <- genes
  
}

cluster_df <- bind_rows(
  lapply(names(cluster_list), function(x){
    
    data.frame(
      Cluster = x,
      Gene = cluster_list[[x]]
    )
    
  })
)

head(cluster_df)

write.csv(cluster_df,file = "human_ratio_cluster_merge.csv",row.names = F)
