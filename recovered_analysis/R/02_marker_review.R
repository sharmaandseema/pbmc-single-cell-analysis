# Recovered marker review; uses the ORIGINAL saved, annotated baseline.
# This does not annotate newly computed clusters.
marker_check <- data.frame(
  Gene = c("TRAC", "TRBC1", "TRBC2", marker_genes))
marker_check$Present <- marker_check$Gene %in% rownames(pbmc_annotated[["RNA"]])
save_table(marker_check, "marker_presence.csv")
missing <- setdiff(marker_genes, rownames(pbmc_annotated[["RNA"]]))
if (length(missing)) stop("Expected plotting markers are missing: ", paste(missing, collapse = ", "))
panels <- list(
  B = c("MS4A1", "CD79A"),
  T = c("CD3D", "CD3E"),
  cytotoxic = c("NKG7", "GNLY", "KLRD1"),
  myeloid = c("LYZ", "CD14", "FCGR3A", "MS4A7"),
  dendritic_platelet = c("FCER1A", "CD1C", "PPBP", "PF4")
)
for (panel_name in names(panels)) {
  marker_violin <- VlnPlot(pbmc_annotated, features = panels[[panel_name]], assay = "RNA",
                          group.by = "seurat_clusters", layer = "data", pt.size = 0.1,
                          add.noise = FALSE, ncol = 2)
  ggplot2::ggsave(file.path(run_dir, "figures", paste0("markers_", panel_name, ".png")),
                  marker_violin, width = 12, height = 7, dpi = 300, bg = "white")
}
qc_groups <- VlnPlot(pbmc_annotated,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  group.by = "seurat_clusters", layer = "counts", pt.size = 0.1,
  add.noise = FALSE, ncol = 3)
ggplot2::ggsave(file.path(run_dir, "figures", "qc_by_original_cluster.png"),
                qc_groups, width = 12, height = 5, dpi = 300, bg = "white")
