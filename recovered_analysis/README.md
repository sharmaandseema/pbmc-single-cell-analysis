# Recovered R source

Run instructions, data provenance and verification status are in the repository's
root README.md. The scripts in this folder retain their recovery provenance.
The numerical pipeline was rerun by the learner on 19 November 2026; the retained
cells and baseline partitions matched exactly. See docs/recovery_provenance.md.

From the project root:

```r
source("recovered_analysis/rebuild_and_check.R")
```

The small historical label-reference CSV supports rebuilding candidate plots
without committing a large Seurat checkpoint. That route stops if cell sets or
partitions differ. The learner executed this route successfully on 19 November
2026 from the prepared GitHub folder; the run completed with matching partitions.
Do not claim that numerical agreement validates the biological annotations.
