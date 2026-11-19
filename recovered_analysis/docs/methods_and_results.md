# PBMC3K: methods, observed results and limitations

**Status, 19 November 2026:** the learner completed the recovered numerical
pipeline after being instructed to restart R. Uploaded Console output confirms
the same retained cell set and exactly matching baseline cluster partitions.
The three sensitivity tables and 51 additional cells also match the earlier
reported results. Source recovery was assistant-supported. The assistant reviewed
the uploaded output and did not execute Seurat locally. The public label-reference
route was subsequently executed successfully from the prepared GitHub folder,
without the historical Seurat RDS, using the existing public count input. The
barcode set and partitions matched before historical candidate labels were attached;
the sensitivity checks and figure exports then completed.

## Aim and public input

Explore PBMC RNA-expression heterogeneity and assess whether broad candidate
annotations are sensitive to selected analysis choices. The public 10x PBMC3K
dataset is from one healthy donor. Input is the Cell Ranger 1.1.0 filtered
gene-barcode matrix, not FASTQ reads and not unfiltered droplets. Upstream cell
calling was already performed by the provider. Rows are genes; columns are
cell barcodes; entries are UMI counts. Cells are not independent donors.

Source: [10x PBMC3K dataset](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0).
[Count archive](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz).
Workflow informed by the [official Seurat PBMC tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).
Provider attribution and applicable data licence should accompany redistribution.

## Historical baseline method

Import yielded **32,738 genes x 2,700 cells**. Seurat object creation retained
all input features (`min.cells = 0`, `min.features = 0`), unlike the tutorial's
initial feature filter. Mitochondrial percentage used counts from gene names
matching `^MT-` divided by each cell's total counts. Baseline retained cells with
**more than 200 and fewer than 2,500 detected genes, and mitochondrial counts
below 5%**. These strict inequalities are dataset-specific teaching choices,
not universal biological quality boundaries. Gene rows were not pruned again.

LogNormalize divided each gene count by cell total counts, multiplied by
10,000 and applied log1p. VST selected 2,000 variable features; scaling and PCA
used those features, with 30 PCs calculated and seed 42. An elbow plot motivated
an exploratory starting choice of 10 PCs. No count-depth or mitochondrial
regression was specified. Neighbors used PCs 1–10 and k = 20. SNN clustering
used Louvain (`algorithm = 1`), resolution 0.5 and random seed 42. UMAP used
PCs 1–10, uwot, 30 neighbors, min.dist 0.3, cosine metric and seed 42 in the
reconstruction. Earlier commands are reconstructed; exact historical default
settings cannot be verified from the truncated upload alone.

## Observed QC and baseline results

| Metric before filtering | Median | Maximum |
|---|---:|---:|
| Total UMI counts | 2,197 | 15,844 |
| Detected genes | 817 | 3,422 |
| Mitochondrial UMI percentage | 2.030% | 22.569% |

| Candidate filter, with the same gene-count bounds | Retained | Excluded |
|---|---:|---:|
| Mitochondrial percentage < 5 | 2,638 | 62 |
| Mitochondrial percentage < 10 | 2,689 | 11 |

There were 2,000 variable features. The reported first ten were PPBP, LYZ,
S100A9, IGLL5, FTL, GNLY, PF4, FTH1, GNG11 and S100A8. They are variable
features, not a ranked list of validated cell-type markers or disease effects.
Nine baseline clusters contained 690, 479, 455, 344, 321, 160, 140, 35 and
14 cells (IDs 0–8). UMAP is an expression-derived visualization; axes are not
physical position, developmental time or evidence of cell interaction.

## Provisional marker-based annotations

| Original cluster(s) | Candidate label | Cells | Selected RNA evidence examined |
|---|---|---:|---|
| 0, 2 | T | 1,145 | CD3D, CD3E |
| 1 | CD14_mono | 479 | LYZ, CD14 |
| 3 | B | 344 | MS4A1, CD79A |
| 4 | Cytotoxic_T | 321 | CD3D/CD3E with cytotoxic-associated NKG7/GNLY |
| 5 | FCGR3A_mono | 160 | LYZ, FCGR3A, MS4A7 |
| 6 | NK_like | 140 | NKG7, GNLY, KLRD1; mixed CD3E signal requires caution |
| 7 | Dendritic_like | 35 | FCER1A, CD1C |
| 8 | Platelet_like | 14 | PPBP, PF4 |

These descriptive counts sum to 2,638 retained cells. They describe the recovered
cell set, not an unbiased estimate of blood population abundance. Canonical
marker associations follow the official Seurat teaching workflow; each label
remains a candidate requiring context and further validation. Cytotoxic RNA
markers do not demonstrate cytotoxic function. The small dendritic-like and
platelet-like groups deserve additional caution. Cluster 6 showed some CD3E
signal: its NK-like label does not establish a homogeneous pure population.
Gene names TRAC, TRBC1 and TRBC2 were absent in the RNA feature list; this does
not establish biological absence of their transcripts. No formal differential
expression tests, reference annotation, doublet calling or ambient correction
were performed. Selected marker plots used to assign labels are not independent
validation of those labels.

## Three sensitivity checks already performed

1. **Resolution 0.5 versus 0.8**, same baseline graph: 9 versus 11 groups.
   Some groups were preserved while T/myeloid-related groups split or mixed.
   This is not a strictly hierarchical subdivision and does not establish
   a uniquely correct resolution. Numeric labels are arbitrary.
2. **10 versus 15 PCs**, k = 20 and resolution 0.5: both produced 9 groups,
   with changes in T-related partitions and some cross-group mixing. For
   example, all 140 original NK-like cells remained together, but that
   alternative group also included 19 original cytotoxic-T candidates.
   B-cell candidates were mostly together (342 of 344 in the main group).
   These comparisons used shared cell barcodes; unchanged cluster count alone
   would not establish stability. The alternative grouping was displayed on
   the existing baseline UMAP, not a newly computed 15-PC UMAP.
3. **Mitochondrial filter < 5% versus < 10%**, same detected-gene bounds and a
   full normalization/HVG/PCA/graph rerun: 51 additional cells were retained.
   Their alternative numeric cluster counts were 15, 11, 3, 5 and 17 in
   clusters 0–4, and zero in 5–8. This is a membership summary, not quality
   validation or biological annotation of the additional cells. The overlap
   table compares the 2,638 shared baseline cells; it excludes the 51 extras.

Broad patterns were retained under several choices, with meaningful boundary
changes. These limited checks support transparent exploratory reporting;
they do not prove the baseline is optimal or every annotation correct.

## Figures and interpretation

The annotated UMAP uses the saved baseline coordinates and manual candidate
labels. The marker dot plot displays normalized RNA, not protein. Dot size
denotes the percentage of cells within each candidate group with detected
expression. Color denotes per-gene scaled group mean expression; color values
cannot be used to compare absolute abundance across different genes. The
grouped average does not show per-cell coexpression. Original plots were
exported as PNG and PDF; the recovery scripts export fresh copies separately.

## Reproducibility record and remaining limits

Learner-reported environment: R 4.6.1, Seurat 6.0.0, SeuratObject 5.4.0,
ggplot2 4.0.4. The learner saved `session_info.txt`; its full contents were not
uploaded. New recovery runs record their own actual software versions.
An annotated RDS checkpoint was confirmed to exist locally. The history is
incomplete. The recovered numerical pipeline has now reproduced the retained
cell set and baseline partition in the learner's reported rerun. Exact UMAP
coordinate equality was not checked. Seeds
alone do not guarantee identical results across package/platform versions.

This is one donor without condition contrasts or independent biological
replication. No donor-level differential expression, disease association,
statistical significance or causal effect is claimed. No batch integration
was needed for a multi-donor contrast because no such contrast was attempted.
Quality filtering can alter the retained composition. Ambient RNA, doublets,
low capture efficiency and uncertain annotations remain possible limitations.
No research novelty, external validation or publication readiness is claimed.
