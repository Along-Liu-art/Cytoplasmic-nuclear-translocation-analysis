# 02_cluster_gene_expression_merge_and_followup.R
#
# Purpose: Cluster gene/expression merging, correlation p-values, and paired cytoplasm/nucleus tables
# Dataset: Human brain nuclear/cytoplasmic proteomics ratio data
# Notes:
#   - Run from the repository root directory.
#   - Input data and intermediate cluster gene lists are not included.
#   - See README.md for expected input files and workflow.
#
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tibble)
})

# ---- 1) 路径与目标 clusters ----
expr_file <- "Human_data_grouped_13ages.csv"
keep_clusters <- c(1, 2, 6, 9, 10)

# ---- 2) 读取表达矩阵（兼容两种格式：Gene为首列 / rownames）----
read_expr <- function(path){
  stopifnot(file.exists(path))
  df1 <- read.csv(path, check.names = FALSE)
  if ("Gene" %in% names(df1)) {
    expr_df <- df1
  } else {
    df2 <- read.csv(path, row.names = 1, check.names = FALSE)
    expr_df <- rownames_to_column(as.data.frame(df2), var = "Gene")
  }
  # 年龄列按数值排序
  age_cols <- setdiff(names(expr_df), "Gene")
  age_cols <- age_cols[order(as.numeric(age_cols))]
  expr_df <- expr_df[, c("Gene", age_cols)]
  expr_df
}

expr_df <- read_expr(expr_file)

# ---- 3) 读取 cluster 基因列表（仅 1/2/6/9/10）----
read_gene_list <- function(fn) {
  if (!file.exists(fn)) stop("找不到文件：", fn)
  lines <- readLines(fn, warn = FALSE)
  # 去空白、去重保序
  uniq <- character(0); seen <- new.env(parent = emptyenv())
  for (s in trimws(lines)) {
    if (nzchar(s) && is.null(seen[[s]])) { uniq <- c(uniq, s); seen[[s]] <- TRUE }
  }
  uniq
}

cluster_files <- setNames(
  paste0("Cluster_", keep_clusters, "_genes.txt"),
  nm = keep_clusters
)

# 若有文件缺失，给提示但继续处理其它
missing_files <- cluster_files[!file.exists(cluster_files)]
if (length(missing_files) > 0) {
  warning("以下 cluster 基因文件缺失，将跳过：",
          paste(names(missing_files), collapse = ", "))
  cluster_files <- cluster_files[file.exists(cluster_files)]
}

gene2cl <- map2_dfr(
  .x = as.integer(names(cluster_files)),
  .y = unname(cluster_files),
  ~ tibble(Gene = read_gene_list(.y), Cluster = .x)
)

if (nrow(gene2cl) == 0) stop("未成功读取任何目标 cluster 的基因。")

# ---- 4) 合并（左连接，保留所有 cluster 基因）----
merged <- gene2cl %>%
  left_join(expr_df, by = "Gene")

# 年龄列顺序（再次确保）
age_cols <- setdiff(names(expr_df), "Gene")
merged  <- merged %>% select(Gene, Cluster, all_of(age_cols))

# ---- 5) 输出：总表 + 每个cluster子表 + 缺失基因 + 汇总 ----
out_combined <- "Merged_Human_13ages_clusters_1_2_6_9_10.csv"
write_csv(merged, out_combined)

# 每个 cluster 的子表（仅保留在矩阵中有任意一个年龄非 NA 的基因）
has_expr <- merged %>%
  mutate(has_expr = rowSums(!is.na(across(all_of(age_cols)))) > 0)

per_cluster_paths <- list()
for (cid in sort(unique(has_expr$Cluster))) {
  sub <- has_expr %>% filter(Cluster == cid, has_expr) %>% select(-has_expr)
  out_path <- sprintf("Merged_cluster_%d_with_expression.csv", cid)
  write_csv(sub, out_path)
  per_cluster_paths[[as.character(cid)]] <- out_path
}

# 缺失基因（所有年龄都是 NA）
missing_tbl <- has_expr %>% filter(!has_expr) %>% select(Gene, Cluster)
if (nrow(missing_tbl) > 0) {
  write_csv(missing_tbl, "Merged_missing_genes_no_expression.csv")
}

# 汇总统计
summary_tbl <- has_expr %>%
  group_by(Cluster) %>%
  summarise(
    n_with_expr = sum(has_expr),
    n_missing   = sum(!has_expr),
    .groups = "drop"
  ) %>% arrange(Cluster)

write_csv(summary_tbl, "Merged_cluster_counts_summary.csv")

# ---- 6) 控制台小结 ----
message("合并完成：")
message(" - 总表：", out_combined)
if (length(per_cluster_paths)) {
  message(" - 子表：")
  walk(names(per_cluster_paths), ~ message("   Cluster ", .x, ": ", per_cluster_paths[[.x]]))
}
if (nrow(missing_tbl) > 0) {
  message(" - 缺失基因列表：Merged_missing_genes_no_expression.csv")
}
message(" - 汇总：Merged_cluster_counts_summary.csv")



#################################计算P值
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

infile  <- "Protein.csv"
outfile <- "Protein_with_pval.csv"

# 读入（保留原列名）
df <- read_csv(infile, show_col_types = FALSE)

# 识别年龄列：列名全是数字的列，如 "21","23",...,"73"
is_age_col <- function(nm) str_detect(nm, "^[0-9]+$")
age_cols <- names(df)[sapply(names(df), is_age_col)]
# 按数值大小排序
age_cols <- age_cols[order(as.numeric(age_cols))]

if (!"Correlation" %in% names(df)) {
  stop("未找到列 'Correlation'，请先在表中提供每行的 Pearson 相关系数。")
}

# 计算每行用于相关的有效样本数（非 NA 的年龄点数）
df <- df %>%
  mutate(
    N_used_for_corr = rowSums(!is.na(across(all_of(age_cols))))
  )

# 根据 r 和 N 计算双侧 p 值
df <- df %>%
  mutate(
    r = as.numeric(Correlation),
    dfree = pmax(N_used_for_corr - 2, 0),                 # 自由度
    Correlation_pval = case_when(
      is.na(r) | N_used_for_corr < 3 ~ NA_real_,          # 样本不足
      abs(r) == 1 & N_used_for_corr >= 3 ~ 0,             # 完美相关
      TRUE ~ {
        tstat <- r * sqrt(dfree / (1 - r^2))
        2 * pt(-abs(tstat), df = dfree)                   # 双侧 p 值
      }
    )
  ) %>%
  select(-r, -dfree)  # 清理临时列

# 将新列追加到末尾并写出
write_csv(df, outfile)
cat("已写出：", outfile, "\n")

# 可选：同时给出 FDR（如需要，取消注释） --------------------
# df$Correlation_fdr <- p.adjust(df$Correlation_pval, method = "BH")
# write_csv(df, outfile)

# 预览
print(df %>% select(Gene, Correlation, N_used_for_corr, Correlation_pval) %>% head())


#######################################合并Cyto和Nuc的数据
# ==========================================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# ---- 1) 读取数据（保留原列名）----
cand <- read_csv("Candidate.csv", show_col_types = FALSE)
cyto <- read_csv("Human_cyto_16.csv", show_col_types = FALSE)
nuc  <- read_csv("Human_Nuc_16.csv",  show_col_types = FALSE)

stopifnot("Gene" %in% names(cand), "Gene" %in% names(cyto), "Gene" %in% names(nuc))

# ---- 2) 构造“配对键”：去掉列名中的 _Cyto / _Nuc 标记（保留复孔编号）----
# 例: 70F_Cyto_1 -> 70F_1,  21F_Cyto -> 21F
mk_base_key <- function(nm) {
  # 把 "_Cyto" 或 "_Nuc" 这一段去掉，并清理多余下划线
  key <- gsub("_(Cyto|Nuc)_?", "_", nm)   # 替换成一个 "_"
  key <- gsub("_+$", "", key)             # 去掉末尾多余 "_"
  key <- gsub("_+", "_", key)             # 压缩连续 "_"
  key
}

cyto_cols <- setdiff(names(cyto), "Gene")
nuc_cols  <- setdiff(names(nuc),  "Gene")

cyto_keys <- tibble(col = cyto_cols, key = mk_base_key(cyto_cols))
nuc_keys  <- tibble(col = nuc_cols,  key = mk_base_key(nuc_cols))

# 只保留两边都存在的配对键（严格配对；含复孔）
paired_keys <- intersect(cyto_keys$key, nuc_keys$key)
if (length(paired_keys) < 3) {
  warning("可配对的样本不足 3 对，相关性将无法稳定估计。")
}

pairs <- tibble(
  key      = paired_keys,
  cyto_col = cyto_keys$col[match(paired_keys, cyto_keys$key)],
  nuc_col  = nuc_keys$col[match(paired_keys,  nuc_keys$key)]
)

# ---- 3) 合并候选基因与 Cyto/Nuc，整理成“宽表”（Gene/Cluster + 成对列）----
cand_basic <- if ("Cluster" %in% names(cand)) cand[, c("Gene","Cluster")] else cand["Gene"]
merged <- cand_basic %>%
  left_join(cyto, by = "Gene") %>%
  left_join(nuc,  by = "Gene")

# 仅保留 Gene/Cluster + 成对列（其余列放到末尾以免打乱）
ordered_cols <- c("Gene", intersect("Cluster", names(merged)))
paired_cols  <- as.vector(t(pairs[, c("cyto_col","nuc_col")]))
ordered_cols <- c(ordered_cols, paired_cols)
leftovers    <- setdiff(names(merged), ordered_cols)
merged_wide  <- merged[, c(ordered_cols, leftovers)]

merged_wide_unique <- merged_wide[!duplicated(merged_wide$Gene), ]

write.csv(merged_wide_unique,file = "merged_wide_unique.csv",row.names = F)
