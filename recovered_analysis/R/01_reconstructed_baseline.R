# RECONSTRUCTED from earlier code and outputs in the conversation.
# Import/QC/PCA/clustering/UMAP are NOT present in the uploaded 512-line history.
# This is a proposed reproducible reconstruction, awaiting a clean-session run.
# It creates numeric clusters; it NEVER assigns cell types by cluster number.
input_dir <- "data/filtered_gene_bc_matrices/hg19"
archive <- "data/pbmc3k_counts.tar.gz"
source_url <- paste0(
  "https://cf.10xgenomics.com/samples/cell/pbmc3k/",
  "pbmc3k_filtered_gene_bc_matrices.tar.gz"
)
if (!dir.exists(input_dir)) {
  if (!file.exists(archive)) {
    stop("Input matrix not found. Keep the original data folder in this project. ",
         "See README.md for the official download command.")
  }
  # Use a temporary extraction location; do not modify the original data.
  extraction_dir <- tempfile("pbmc3k_extract_")
  dir.create(extraction_dir)
  untar(archive, exdir = extraction_dir)
  input_dir <- file.path(extraction_dir, "filtered_gene_bc_matrices", "hg19")
  if (!dir.exists(input_dir)) stop("The archive has an unexpected folder layout.")
}
input_files <- list.files(input_dir, full.names = TRUE)
writeLines(c("Official source:", source_url, "Files used:", input_files),
           file.path(run_dir, "logs", "input_source.txt"))
checksums <- tools::md5sum(input_files)
save_table(data.frame(File = basename(names(checksums)), MD5 = unname(checksums)),
           "input_checksums.csv")
# Checksums document the local files; no provider-authenticated checksum is claimed.
pbmc_counts <- Read10X(data.dir = input_dir)
if (!is.matrix(pbmc_counts) && !inherits(pbmc_counts, "Matrix")) {
  stop("Expected one gene-expression matrix from this PBMC3K archive.")
}
print(dim(pbmc_counts))
pbmc <- CreateSeuratObject(
  counts = pbmc_counts, project = "PBMC3K", min.cells = 0, min.features = 0
)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
before_meta <- pbmc[[]]
pass5 <- with(before_meta, nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pass10 <- with(before_meta, nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
filter_comparison <- data.frame(
  Candidate = c("Mito below 5%", "Mito below 10%"),
  Retained = c(sum(pass5), sum(pass10)),
  Excluded = c(sum(!pass5), sum(!pass10))
)
print(filter_comparison)
save_table(filter_comparison, "qc_filter_comparison.csv")
save_table(data.frame(Cell = rownames(before_meta), before_meta, check.names = FALSE),
           "qc_metrics_before_filtering.csv")
qc_plot <- VlnPlot(pbmc, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                  pt.size = 0.1, add.noise = FALSE, ncol = 3)
ggplot2::ggsave(file.path(run_dir, "figures", "qc_before_filtering.png"),
                qc_plot, width = 12, height = 5, dpi = 300, bg = "white")
pbmc_qc <- subset(pbmc, cells = colnames(pbmc)[pass5])
pbmc_norm <- NormalizeData(pbmc_qc, assay = "RNA",
                          normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_norm <- FindVariableFeatures(pbmc_norm, assay = "RNA",
                                 selection.method = "vst", nfeatures = 2000)
variable_genes <- VariableFeatures(pbmc_norm)
save_table(data.frame(Gene = variable_genes), "variable_genes.csv")
pbmc_pca <- ScaleData(pbmc_norm, assay = "RNA", features = variable_genes)
pbmc_pca <- RunPCA(pbmc_pca, assay = "RNA", features = variable_genes,
                   npcs = 30, seed.use = 42)
elbow <- ElbowPlot(pbmc_pca, ndims = 30)
ggplot2::ggsave(file.path(run_dir, "figures", "pca_elbow.png"),
                elbow, width = 8, height = 5, dpi = 300, bg = "white")
set.seed(42)
pbmc_clustered <- FindNeighbors(pbmc_pca, reduction = "pca", dims = 1:10, k.param = 20)
pbmc_clustered <- FindClusters(pbmc_clustered, graph.name = "RNA_snn",
                               resolution = 0.5, algorithm = 1, random.seed = 42)
pbmc_umap <- RunUMAP(pbmc_clustered, reduction = "pca", dims = 1:10,
                     umap.method = "uwot", n.neighbors = 30, min.dist = 0.3,
                     metric = "cosine", seed.use = 42)
saveRDS(pbmc_umap, file.path(run_dir, "pbmc_rebuilt_numeric_clusters.rds"))
save_table(data.frame(Cell = colnames(pbmc_umap),
                      Rebuilt_cluster = as.character(pbmc_umap$seurat_clusters)),
           "rebuilt_cell_clusters.csv")
numeric_umap <- DimPlot(pbmc_umap, reduction = "umap", group.by = "seurat_clusters",
                        label = TRUE) + ggplot2::labs(title = "Rebuilt numeric clusters")
ggplot2::ggsave(file.path(run_dir, "figures", "rebuilt_numeric_umap.png"),
                numeric_umap, width = 10, height = 7, dpi = 300, bg = "white")
if (file.exists(checkpoint)) {
  historical <- load_historical_baseline()
  same_cells <- setequal(colnames(historical), colnames(pbmc_umap))
  check_lines <- paste("Same retained cell set:", same_cells)
  if (same_cells) {
    # Align by barcode, not current row order. Cluster numbers can be permuted.
    aligned <- pbmc_umap[[]][colnames(historical), "seurat_clusters"]
    historical_vs_rebuilt <- table(
      Historical = historical$seurat_clusters, Rebuilt = aligned)
    save_table(as.data.frame(historical_vs_rebuilt), "historical_vs_rebuilt.csv")
    same_partition <- all(rowSums(historical_vs_rebuilt > 0) == 1) &&
                      all(colSums(historical_vs_rebuilt > 0) == 1)
    check_lines <- c(check_lines, paste("Exactly matching cell partitions:", same_partition))
  }
  writeLines(c(check_lines,
    "Cell-type labels were not automatically assigned to rebuilt clusters.",
    "Review marker evidence before interpreting any changed groups."),
    file.path(run_dir, "logs", "baseline_rebuild_check.txt"))
  message(paste(check_lines, collapse = "\n"))
}
