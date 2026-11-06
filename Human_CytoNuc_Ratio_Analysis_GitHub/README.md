# Human Brain Cytoplasm–Nucleus Proteomics Ratio Analysis

R scripts for analyzing age-associated changes in the nuclear/cytoplasmic protein ratio in human brain samples.

> **Scope:** Human ratio data only. Mouse analyses and raw/proprietary proteomics data are not included.

## Repository structure

```text
.
├── README.md
├── .gitignore
├── data/
│   └── input/       # Place input CSVs here (not tracked by default)
├── scripts/
│   ├── 01_human_ratio_trajectory_clustering.R
│   ├── 02_cluster_gene_expression_merge_and_followup.R
│   ├── 03_cluster_amplitude_variability.R
│   ├── 04_kegg_enrichment_and_heatmap.R
│   ├── 05_export_cluster_gene_table.R
│   └── 06_cluster_tables_and_old_young_ratios.R
└── results/          # Generated files (not tracked by default)
```

## Analysis workflow

1. **Ratio trajectory clustering** — `01_human_ratio_trajectory_clustering.R`
   - Reads age-grouped human ratio matrix.
   - Orders age columns numerically, linearly interpolates missing values where possible, removes proteins that remain incomplete or have zero variance, and row-wise Z-scores the data.
   - Performs Euclidean-distance hierarchical clustering with Ward.D2 and cuts the tree into 12 clusters.
   - Exports cluster gene lists, dendrograms, and cluster-average LOESS trajectories.
2. **Cluster table / follow-up data preparation** — `02_cluster_gene_expression_merge_and_followup.R`
   - Merges selected cluster gene lists with expression tables.
   - Calculates two-sided correlation p-values from a supplied Pearson correlation and the number of non-missing age points.
   - Pairs cytoplasmic and nuclear sample columns using sample-name conventions and exports a wide table.
3. **Cluster amplitude and variability** — `03_cluster_amplitude_variability.R`
   - Computes per-cluster mean protein amplitude (per-protein max–min, then averaged) and mean within-cluster SD across ages.
4. **KEGG enrichment and heatmap** — `04_kegg_enrichment_and_heatmap.R`
   - Converts gene symbols to Entrez IDs, performs KEGG enrichment, and exports cluster-specific and combined results.
5. **Cluster membership export** — `05_export_cluster_gene_table.R`
   - Writes cluster membership in long format.
6. **Cluster tables and old/young ratios** — `06_cluster_tables_and_old_young_ratios.R`
   - Exports combined cluster membership and calculates old/young ratios for selected clusters.






## Citation

If these scripts support a publication, cite the associated study and provide the repository DOI/version when available.
