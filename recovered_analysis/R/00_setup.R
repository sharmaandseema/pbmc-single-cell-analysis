# Shared setup. Run using one of the two runner scripts in this folder.
# Seurat is already installed on the learner's computer. If it is missing:
# install.packages("Seurat")  # Console, once; wait for the > prompt.
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop('Seurat is missing. In the Console run install.packages("Seurat").')
}
library(Seurat)
if (!file.exists("PBMC_project.Rproj")) {
  stop("Open PBMC_project.Rproj first, then run this script from that project.")
}
# A NEW output folder for every run; existing analysis outputs stay intact.
dir.create("results", showWarnings = FALSE)
run_dir <- tempfile(
  pattern = paste0("recovered_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
  tmpdir = "results"
)
dir.create(run_dir)
for (folder in c("figures", "tables", "logs")) {
  dir.create(file.path(run_dir, folder))
}
save_table <- function(x, name) {
  write.csv(x, file.path(run_dir, "tables", name), row.names = FALSE)
}
checkpoint <- "data/processed/pbmc_annotated_candidates.rds"
marker_genes <- c(
  "CD3D", "CD3E", "MS4A1", "CD79A", "NKG7", "GNLY", "KLRD1",
  "LYZ", "CD14", "FCGR3A", "MS4A7", "FCER1A", "CD1C", "PPBP", "PF4"
)
load_historical_baseline <- function() {
  if (file.exists(checkpoint)) {
    object <- readRDS(checkpoint)
  } else {
    reference_file <- "reference/historical_cell_labels.csv"
    if (!file.exists(reference_file) || !exists("pbmc_umap", inherits = TRUE)) {
      stop("Run rebuild_and_check.R from the project root to rebuild counts and use the label reference.")
    }
    reference <- read.csv(reference_file, stringsAsFactors = FALSE)
    needed <- c("Cell", "Original_cluster", "Candidate")
    if (!all(needed %in% names(reference)) || anyNA(reference[needed]) ||
        anyDuplicated(reference$Cell) ||
        !setequal(reference$Cell, colnames(pbmc_umap))) {
      stop("Reference labels and rebuilt cells do not match; review input and QC.")
    }
    reference <- reference[match(colnames(pbmc_umap), reference$Cell), ]
    agreement <- table(Reference = reference$Original_cluster,
                       Rebuilt = pbmc_umap$seurat_clusters)
    same_partition <- all(rowSums(agreement > 0) == 1) &&
                      all(colSums(agreement > 0) == 1)
    save_table(as.data.frame(agreement), "reference_vs_rebuilt.csv")
    if (!same_partition) {
      stop("Rebuilt partitions differ from the reference. No candidate labels were transferred. Review markers.")
    }
    object <- pbmc_umap
    object$celltype_candidate <- reference$Candidate
    # Preserve historical IDs as display metadata after checking exact partitions.
    object$seurat_clusters <- as.character(reference$Original_cluster)
    Idents(object) <- "seurat_clusters"
    writeLines(c("Candidate labels copied by barcode from historical marker-based review.",
                 "Exact cell-set and cluster-partition agreement checked before transfer.",
                 "These labels are not independent biological validation."),
               file.path(run_dir, "logs", "annotation_reference_check.txt"))
    message("Reference cell set and cluster partitions match; historical candidate labels attached by barcode.")
  }
  if (!inherits(object, "Seurat") ||
      !all(c("seurat_clusters", "celltype_candidate") %in% colnames(object[[]])) ||
      !"umap" %in% Reductions(object)) {
    stop("The saved object does not contain the expected annotations and UMAP.")
  }
  if (anyNA(object$celltype_candidate)) stop("Some saved candidate labels are missing.")
  DefaultAssay(object) <- "RNA"
  object
}
