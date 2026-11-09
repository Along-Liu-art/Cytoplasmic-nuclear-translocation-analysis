# 01_human_ratio_trajectory_clustering.R
setwd("E:/Cai lab data/Proteomics data/Renew_2025_07/Human/处理过数据")
######加载包
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(tidyr)
library(tibble)
library(tidyverse)
#######读取数据
dat <- read.csv(
  "Human_Ratio_16.csv",
  check.names = FALSE
)

gene <- dat$Gene
age <- as.numeric(colnames(dat)[-1])

ratio <- as.matrix(dat[, -1])
storage.mode(ratio) <- "numeric"

# QC
n_nonmissing <- rowSums(!is.na(ratio))
row_sd <- apply(ratio, 1, sd, na.rm = TRUE)

keep <- n_nonmissing >= 11 &
  is.finite(row_sd) &
  row_sd > 0

ratio <- ratio[keep, , drop = FALSE]
gene <- gene[keep]

# -------------------------
# Z-score
# -------------------------
z_mat <- t(apply(ratio, 1, function(x) {
  (x - mean(x, na.rm = TRUE)) /
    sd(x, na.rm = TRUE)
}))

# -------------------------
# LOESS
# -------------------------
trajectory_mat <- matrix(
  NA,
  nrow = nrow(z_mat),
  ncol = length(age)
)

trajectory_mat <- matrix(
  NA_real_,
  nrow = nrow(z_mat),
  ncol = length(age)
)

rownames(trajectory_mat) <- rownames(z_mat)

for (i in seq_len(nrow(z_mat))) {
  
  y <- z_mat[i, ]
  
  ok <- is.finite(y) & is.finite(age)
  
  # 至少需要5个有效点
  if (sum(ok) < 5) {
    next
  }
  
  dat_loess <- data.frame(
    age = age[ok],
    z = y[ok]
  )
  
  fit <- tryCatch(
    loess(
      z ~ age,
      data = dat_loess,
      span = 0.75,##取大于75%的数据点
      degree = 2,#取二次函数计算
      family = "gaussian",#####使用差值平均方法
      control = loess.control(surface = "direct")
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    next
  }
  
  pred <- predict(
    fit,
    newdata = data.frame(age = age)
  )
  
  trajectory_mat[i, ] <- pred
}


keep_trajectory <- apply(
  trajectory_mat,
  1,
  function(x) all(is.finite(x))
)

table(keep_trajectory)

trajectory_mat <- trajectory_mat[
  keep_trajectory,
  ,
  drop = FALSE
]


rownames(trajectory_mat) <- gene

################################################################分析cluster分离度################
library(cluster)

sil <- silhouette(
  cluster,
  dist(z_mat)
)

str(sil)
mean(sil[, 3])


d <- dist(z_mat)

k_values <- 2:30

sil_results <- data.frame(
  k = k_values,
  mean_silhouette = NA_real_
)

for (i in seq_along(k_values)) {
  
  k <- k_values[i]
  
  cl <- cutree(hc, k = k)
  
  sil <- silhouette(cl, d)
  
  sil_results$mean_silhouette[i] <- mean(sil[, "sil_width"])
}

print(sil_results)
#############################################
# -------------------------
# Euclidean
# -------------------------
d <- dist(
  trajectory_mat,
  method = "euclidean"
)

# -------------------------
# Ward.D2
# -------------------------
hc <- hclust(
  d,
  method = "ward.D2"
)

# -------------------------
# 13 clusters
# -------------------------
cluster <- cutree(
  hc,
  k = 13
)
cluster_info <- data.frame(
  Gene = rownames(trajectory_mat),
  Cluster = factor(cluster)
)
head(cluster_info)
table(cluster_info$Cluster)


#######准备绘图
colnames(trajectory_mat) <- colnames(z_mat)

plot_df <- as.data.frame(trajectory_mat)

plot_df$Gene <- rownames(plot_df)
plot_df$Cluster <- factor(cluster)

plot_df <- plot_df %>%
  pivot_longer(
    cols = -c(Gene, Cluster),
    names_to = "Age",
    values_to = "Z"
  ) %>%
  mutate(
    Age = as.numeric(Age)
  )

print(plot_df)

cluster_summary <- plot_df %>%
  group_by(Cluster, Age) %>%
  summarise(
    mean_Z = mean(Z, na.rm = TRUE),
    median_Z = median(Z, na.rm = TRUE),
    q10 = quantile(Z, 0.10, na.rm = TRUE),
    q90 = quantile(Z, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

cluster_cor <- cluster_summary %>%
  group_by(Cluster) %>%
  summarise(
    
    rho = cor(
      Age,
      mean_Z,
      method = "spearman"
    ),
    
    p = cor.test(
      Age,
      mean_Z,
      method = "spearman",
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  
  mutate(
    FDR = p.adjust(p, method = "BH")
  )

print(cluster_cor)
cluster_cor <- cluster_cor %>%
  mutate(
    Cluster = as.character(Cluster),
    Cluster_num = as.numeric(Cluster)
  )

# cluster size
cluster_n <- cluster_info %>%
  count(Cluster, name = "N")
cluster_n <- cluster_n %>%
  mutate(
    Cluster = as.character(Cluster)
  )
head(cluster_n)
# 合并 cluster size
cluster_cor <- cluster_cor %>%
  left_join(
    cluster_n,
    by = "Cluster"
  ) %>%
  mutate(
    facet_label = paste0(
      "Cluster ", Cluster,
      "\nn = ", N
    ),
    facet_label = factor(
      facet_label,
      levels = facet_label[order(Cluster_num)]
    )
  )


cluster_summary <- cluster_summary %>%
  left_join(
    cluster_cor %>%
      dplyr::select(
        Cluster,
        facet_label
      ),
    by = "Cluster"
  )


cluster_summary <- cluster_summary %>%
  mutate(
    Cluster = as.character(Cluster),
    Cluster_num = as.numeric(Cluster)
  )


head(cluster_summary)

write.csv(cluster_summary,file = "cluster_summary.csv",row.names = F)
write.csv(cluster_cor,file = "cluster_cor.csv",row.names = F)


######################################################################################
pdf("Humanratio_Protein_Trajectories_by_Cluster.pdf",height = 7,width = 7)
p <- ggplot(
  cluster_summary,
  aes(x = Age, y = mean_Z)
) +
  
  # cluster内部 variation
  geom_ribbon(
    aes(
      ymin = q10,
      ymax = q90
    ),
    fill = "grey75",
    alpha = 0.6
  ) +
  
  # cluster平均trajectory
  geom_line(
    linewidth = 1.2,
    color = "black"
  ) +
  
  # zero line
  #geom_hline(
   # yintercept = 0,
    #linewidth = 0.3,
    #color = "grey60"
  #) +
  facet_wrap(
    ~ facet_label,
    ncol = 4
  ) +
  
  scale_x_continuous(
    breaks = c(20, 40, 60, 80)
  ) +
  
  coord_cartesian(
    ylim = c(-2, 2)
  ) +
  
  labs(
    x = "Age",
    y = "Z-score"
  ) +
  
  theme_classic() +
  
  theme(
    strip.background = element_rect(
      fill = "grey90",
      color = "black"
    ),
    strip.text = element_text(
      color = "black"
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    panel.grid = element_blank()
  )

p

dev.off()

rownames(z_mat) <- gene
head(trajectory_mat,n=6)
head(z_mat,n=6)


####################################################################
###################################################################################################################################################

#########聚类
head(cluster_info)
table(cluster_info$Cluster)
gene_map <- bitr(
  cluster_info$Gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

cluster_gene <- cluster_info %>%
  inner_join(gene_map, by = c("Gene" = "SYMBOL")) %>%
  distinct(Cluster, ENTREZID)

clusters <- sort(unique(cluster_gene$Cluster))

kegg_list <- list()

for (cl in clusters) {
  
  genes <- cluster_gene %>%
    filter(Cluster == cl) %>%
    pull(ENTREZID) %>%
    unique()
  
  ek <- enrichKEGG(
    gene = genes,
    organism = "hsa",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500
  )
  
  if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
    
    tmp <- as.data.frame(ek)
    tmp$Cluster <- cl
    
    kegg_list[[as.character(cl)]] <- tmp
  }
}

kegg_all <- bind_rows(kegg_list)

head(kegg_all)


kegg_top5 <- kegg_all %>%
  filter(p.adjust < 0.05) %>%
  group_by(Cluster) %>%
  arrange(p.adjust, pvalue, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()


kegg_top5 %>%
  dplyr::select(Cluster, Description, p.adjust, Count)

plot_kegg <- kegg_top5 %>%
  mutate(
    neg_log10_FDR = -log10(p.adjust)
  )

plot_kegg %>%
  dplyr::select(Cluster, Description, p.adjust, neg_log10_FDR)

pathway_order <- plot_kegg %>%
  group_by(Description) %>%
  summarise(
    best_FDR = min(p.adjust),
    .groups = "drop"
  ) %>%
  arrange(best_FDR) %>%
  pull(Description)

plot_kegg <- plot_kegg %>%
  mutate(
    Description = factor(
      Description,
      levels = rev(pathway_order)
    ),
    Cluster = factor(
      Cluster,
      levels = clusters
    )
  )
plot_kegg <- na.omit(plot_kegg)
# === Step 6: 生成矩阵并补齐缺失 cluster ==

# 统一 cluster 名称
plot_kegg2 <- plot_kegg %>%
  mutate(
    Cluster = paste0("Cluster_", Cluster)
  )


mat_df <- plot_kegg2 %>%
  dplyr::select(Description, Cluster, neg_log10_FDR) %>%
  pivot_wider(
    names_from = Cluster,
    values_from = neg_log10_FDR,
    values_fill = 0
  )



# 补齐缺失的 cluster 列
desired_order <- paste0("Cluster_", 1:13)

missing_clusters <- setdiff(
  desired_order,
  colnames(mat_df)
)

if (length(missing_clusters) > 0) {
  mat_df[missing_clusters] <- 0
}

mat_df <- mat_df %>%
  dplyr::select(
    Description,
    all_of(desired_order)
  )

mat <- mat_df %>%
  column_to_rownames("Description") %>%
  as.matrix()
# === Step 7: 标记关键通路 ===

highlight_terms <- c("Nucleocytoplasmic transport", "Spliceosome", "mRNA surveillance pathway")
highlight_rows <- rownames(mat) %>%
  str_to_lower() %>%
  str_detect(paste(tolower(highlight_terms), collapse = "|"))

# === Step 8: 绘图 ===
library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(c(0, max(mat)), c("#f2f2f2", "#4575b4"))
pdf("KEGG_Top5_Heatmap_R_fixed_20260930_humanratio.pdf", width = 7, height = 7)
ht <- Heatmap(
  mat,
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
draw(ht)
dev.off()


kegg_list <- list()

for (cl in clusters) {
  
  cat("Processing Cluster", cl, "...\n")
  
  genes <- cluster_gene %>%
    filter(Cluster == cl) %>%
    pull(ENTREZID) %>%
    unique()
  
  if (length(genes) < 10) {
    cat("  Too few genes, skipped.\n")
    next
  }
  
  # ---------------------------------
  # Retry up to 3 times
  # ---------------------------------
  
  ek <- NULL
  
  for (attempt in 1:3) {
    
    cat("  Attempt", attempt, "...\n")
    
    ek <- tryCatch(
      
      enrichKEGG(
        gene = genes,
        organism = "hsa",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 10,
        maxGSSize = 500
      ),
      
      error = function(e) {
        cat(
          "  Error:",
          conditionMessage(e),
          "\n"
        )
        NULL
      }
    )
    
    if (!is.null(ek)) {
      break
    }
    
    Sys.sleep(5)
  }

##################################################################
