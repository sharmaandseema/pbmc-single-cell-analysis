# Cleaned from the uploaded history. Original saved baseline is the reference.
# Requires the 2,700-cell pbmc object from 01_reconstructed_baseline.R for mt10.
# The resolution and PC alternatives use the SAME original UMAP coordinates
# for display; no new UMAP is computed for those alternatives.
pbmc_res08 <- FindClusters(pbmc_annotated, graph.name = "RNA_snn",
                           resolution = 0.8, algorithm = 1, random.seed = 42)
baseline_cells <- colnames(pbmc_annotated)
resolution_comparison <- table(
  Original_0.5 = pbmc_annotated$seurat_clusters,
  Alternative_0.8 = pbmc_res08[[]][baseline_cells, "seurat_clusters"])
save_table(as.data.frame(resolution_comparison), "resolution_comparison.csv")
set.seed(42)
pbmc_pc15 <- FindNeighbors(pbmc_annotated, reduction = "pca", dims = 1:15, k.param = 20)
pbmc_pc15 <- FindClusters(pbmc_pc15, graph.name = "RNA_snn",
                          resolution = 0.5, algorithm = 1, random.seed = 42)
pc_comparison <- table(
  Original_10PC = pbmc_annotated$seurat_clusters,
  Alternative_15PC = pbmc_pc15[[]][baseline_cells, "seurat_clusters"])
save_table(as.data.frame(pc_comparison), "pc_comparison.csv")
pbmc_mt10 <- subset(pbmc,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
pbmc_mt10 <- NormalizeData(pbmc_mt10,
  normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_mt10 <- FindVariableFeatures(pbmc_mt10, selection.method = "vst", nfeatures = 2000)
pbmc_mt10 <- ScaleData(pbmc_mt10, features = VariableFeatures(pbmc_mt10))
pbmc_mt10 <- RunPCA(pbmc_mt10, features = VariableFeatures(pbmc_mt10), npcs = 30, seed.use = 42)
set.seed(42)
pbmc_mt10 <- FindNeighbors(pbmc_mt10, reduction = "pca", dims = 1:10, k.param = 20)
pbmc_mt10 <- FindClusters(pbmc_mt10, graph.name = "RNA_snn",
                          resolution = 0.5, algorithm = 1, random.seed = 42)
if (!all(baseline_cells %in% colnames(pbmc_mt10))) {
  stop("The mt10 run does not include every historical baseline cell; inspect input/QC.")
}
mt_comparison <- table(
  Baseline_5pct = pbmc_annotated$seurat_clusters,
  Alternative_10pct = pbmc_mt10[[]][baseline_cells, "seurat_clusters"])
save_table(as.data.frame(mt_comparison), "mitochondrial_filter_comparison.csv")
added_cells <- setdiff(colnames(pbmc_mt10), baseline_cells)
added_cluster_counts <- table(pbmc_mt10[[]][added_cells, "seurat_clusters"])
save_table(data.frame(Cluster = names(added_cluster_counts),
                      Cells = as.integer(added_cluster_counts)), "additional_cell_clusters.csv")
save_table(data.frame(Cell = added_cells,
  Alternative_cluster = as.character(pbmc_mt10[[]][added_cells, "seurat_clusters"])),
  "additional_cell_assignments.csv")
print(resolution_comparison)
print(pc_comparison)
print(mt_comparison)
print(length(added_cells))
for (label in c("resolution08", "pc15")) {
  object <- if (label == "resolution08") pbmc_res08 else pbmc_pc15
  alternative_plot <- DimPlot(object, reduction = "umap", group.by = "seurat_clusters",
                               label = TRUE, repel = TRUE) +
    ggplot2::labs(title = paste("Alternative grouping:", label),
                  caption = "Shown on the original baseline UMAP; coordinates were not recomputed.")
  ggplot2::ggsave(file.path(run_dir, "figures", paste0(label, "_original_umap.png")),
                  alternative_plot, width = 11, height = 7, dpi = 300, bg = "white")
}
