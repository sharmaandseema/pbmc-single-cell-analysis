# Recovery provenance

- Uploaded evidence: `pbmc_console_history.txt`, 512 lines.
- SHA-256 of uploaded bytes: `86354fb42f955bd6fdca3435852128449e21bc1926733d254ad2f9615de16733`.
- First recorded command: cytotoxic-marker `VlnPlot` of NKG7/GNLY.
- History documents later commands, not their Console outputs. Historical
  numerical results in the report come from the learner's messages/screenshots.
- Initial import through baseline UMAP: reconstructed from conversation;
  absent from the upload. Exact original commands/defaults remain unverified.
- Later marker review, annotation, sensitivity checks, figures and exports:
  cleaned from upload, with barcode alignment, explicit dependencies, protected
  paths and separate outputs added for recovery.
- Added missing dependency: define `marker_genes` before using DotPlot.
- Removed repeated attempts, absent TRAC plotting request and stray `width = 14,`.
- No attempt was made to execute the uploaded history as a script.
- No access to the learner's local computer or raw matrix/checkpoint was used.
- R is not installed in the assistant workspace; these recovered R scripts
  have not been run in R here. R syntax was checked with a tree-sitter parser,
  and official Seurat interfaces were checked. Clean-session R execution remains
  the next validation step.
- Source recovery was assistant-supported; project learning and interpretation
  should be described transparently in a future portfolio README.

Official documentation consulted 19 November 2026:
- https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
- https://satijalab.org/seurat/reference/findclusters
- https://satijalab.org/seurat/reference/dotplot

Historical manual mapping, recorded here for audit only:
0=T; 1=CD14_mono; 2=T; 3=B; 4=Cytotoxic_T; 5=FCGR3A_mono;
6=NK_like; 7=Dendritic_like; 8=Platelet_like.
This mapping is NOT applied automatically to reconstructed cluster numbers.

## Subsequent learner verification, 19 November 2026

The uploaded rebuild Console output reports:

- Input dimensions: 32,738 genes by 2,700 cells.
- Baseline retained: 2,638 cells; alternative mt filter retained: 2,689.
- Same retained cell set: TRUE.
- Exactly matching cell partitions: TRUE.
- Baseline: 9 groups; resolution 0.8 alternative: 11 groups.
- PC and mitochondrial comparison tables match earlier learner-reported tables.
- Additional cells under the 10% mitochondrial filter: 51.
- Final completion message and unchanged candidate counts are present.

The earlier statement that clean-session execution remains outstanding is
superseded by this verification for the recovered checkpoint-based runner.
The user was instructed to restart R first; the uploaded output does not itself
prove the absence of every prior in-memory object. The local runner uses its own
scope and recreates its numerical inputs. Exact UMAP coordinate equality was not
tested. Historical labels remain provisional and were not independently validated.

The GitHub package adds a new reference-based annotation route for clones without
the historical RDS. It requires exactly the same cell set and one-to-one matching
cluster partitions before copying historical labels by barcode. Otherwise it stops.
The learner subsequently executed this route successfully on 19 November 2026.
The uploaded output contains the reference-match message, all three unchanged
sensitivity tables, 51 additional cells, unchanged candidate counts and a completion
message. The run was from the github_ready folder without the historical RDS,
using a local copy of the public count input. This does not test another machine,
a fresh data download or cross-version reproducibility. Numerical and sensitivity
algorithms were preserved.
