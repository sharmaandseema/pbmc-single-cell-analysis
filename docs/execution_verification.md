# Execution verification: 19 November 2026

Two learner-executed runs were reviewed from uploaded Console output:

1. Recovered source, with the historical saved baseline available: same retained
   cell set TRUE and exactly matching cell partitions TRUE; all sensitivity checks
   and final exports completed.
2. Prepared GitHub project, without the historical RDS: reference cell set and
   cluster partitions matched; historical candidate labels attached by barcode;
   sensitivity checks and final exports completed.

Input dimensions were 32,738 genes by 2,700 cells. The baseline retained 2,638
cells; the alternative mitochondrial filter retained 2,689, adding 51 cells.
Baseline clustering produced 9 groups and the higher-resolution check 11 groups.
Candidate counts and all three contingency tables matched earlier reported results.

These checks demonstrate local execution and numerical partition agreement.
They do not independently validate the candidate cell-type annotations, establish
identical UMAP coordinates, or establish reproducibility across machines/versions.
The public-label route reused a local copy of the public input; it did not test
a fresh network download. No formal biological significance testing was performed.

Observed warnings: Seurat, SeuratObject and sp were built under R 4.6.2 while the
reported session used R 4.6.1; feature underscores were replaced with dashes;
the pre-normalization QC plotting call used counts because no data layer existed.
The observed runs nevertheless completed and passed the partition comparison.

The committed figures and tables record the earlier verified run; both checks
reproduced the reported candidate counts and comparison tables. New run outputs
remain in ignored results/recovered_* folders. Public labels remain provisional.
