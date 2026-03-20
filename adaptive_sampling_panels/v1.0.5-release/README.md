# v1.0.5-20260316-release

This folder contains the `v1.0.5-20260316-release` adaptive sampling release generated from the versioned three-step workflow:

1. normalize inputs to a minimal BED
2. generate a padded panel BED/XLSX
3. merge nearby intervals into a final panel BED/XLSX

## Files
- `v1.0.5-normalized-inputs.bed`
- `v1.0.5-pre-merge-panel-20kb-pad.bed`
- `v1.0.5-pre-merge-panel-20kb-pad.xlsx`
- `v1.0.5-sequencing-panel-20kb-pad-50kb-merge.bed`
- `v1.0.5-sequencing-panel-20kb-pad-50kb-merge.xlsx`
- `v1.0.5-analysis-set-0kb-pad.bed.bed`
- `v1.0.5-analysis-set-0kb-pad.bed.xlsx`

## Difference from previous version
Compared with `v1.0.4-20251224`, this release expands the target set while keeping the same three manual loci (`IGH`, `TRB`, `TRAV`):
- adds the extended-panel genes:
  `B2M`, `CIITA`, `IL4R`, `MYD88`, `PDCD1LG2`, `RCC1`, `RIC1`, `RMI2`, `SNHG3`, `SOCS1`
- adds `CEBPE`
- adds `TFE3`
- the resulting source gene list contains 392 genes plus 3 loci

## Workflow notes
- The normalized BED represents the inputs normalized to coordinate form.
- The analysis set represents the coordinate-ordered regions without padding for downstream targeted analysis.
- The pre-merge panel represents the coordinate-ordered regions with the designated padding.
- The sequencing-panel is the pre-merge panel with regions within the threshold merged. This is the file that should be served to MinKnow for adaptive sampling sequencing runs.
