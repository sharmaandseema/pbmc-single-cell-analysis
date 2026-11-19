# Uses historical candidate labels: checkpoint or partition-checked reference.
candidate_counts <- as.data.frame(table(pbmc_annotated$celltype_candidate))
colnames(candidate_counts) <- c("Candidate", "Cells")
save_table(candidate_counts, "candidate_cell_counts.csv")
save_table(data.frame(Cell = colnames(pbmc_annotated),
  Original_cluster = as.character(pbmc_annotated$seurat_clusters),
  Candidate = pbmc_annotated$celltype_candidate), "historical_cell_labels.csv")
writeLines(c("Figure source: historical checkpoint, or rebuilt baseline with checked historical labels",
             if (file.exists(checkpoint)) checkpoint else "reference/historical_cell_labels.csv",
             "These are historical manual candidate labels, not newly validated identities."),
           file.path(run_dir, "logs", "figure_provenance.txt"))
umap_final <- DimPlot(pbmc_annotated, reduction = "umap", group.by = "celltype_candidate",
  label = TRUE, repel = TRUE, pt.size = 0.3) +
  ggplot2::labs(title = "PBMC3K: provisional cell annotations",
    subtitle = paste(format(ncol(pbmc_annotated), big.mark = ","),
                     "retained cells | historical baseline: 10 PCs; resolution 0.5"),
    caption = "One donor; marker-based candidate labels.") +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 16),
                 plot.subtitle = ggplot2::element_text(size = 11))
missing <- setdiff(marker_genes, rownames(pbmc_annotated[["RNA"]]))
if (length(missing)) stop("Missing marker genes: ", paste(missing, collapse = ", "))
marker_plot <- DotPlot(pbmc_annotated, features = marker_genes, assay = "RNA",
  group.by = "celltype_candidate", cols = c("lightgrey", "blue"),
  dot.scale = 6, scale = TRUE) + RotatedAxis() +
  ggplot2::labs(title = "PBMC3K: marker RNA evidence", x = "Gene", y = "Candidate group",
                caption = "Annotation markers; not independent validation.") +
  ggplot2::guides(
    colour = ggplot2::guide_colourbar(title = "Scaled mean\nexpression"),
    size = ggplot2::guide_legend(title = "Cells detected (%)")) +
  ggplot2::theme(plot.title = ggplot2::element_text(size = 14),
    plot.caption = ggplot2::element_text(size = 9),
    axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 10),
    legend.title = ggplot2::element_text(size = 10))
for (extension in c("png", "pdf")) {
  ggplot2::ggsave(file.path(run_dir, "figures", paste0("umap_candidates.", extension)),
                  umap_final, width = 12, height = 7, units = "in", dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(run_dir, "figures", paste0("marker_evidence.", extension)),
                  marker_plot, width = 14, height = 8, units = "in", dpi = 300, bg = "white")
}
print(candidate_counts)
