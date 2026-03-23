# Gene Panels

Repository for versioned gene panel inputs and adaptive sampling panel outputs used in sequencing workflows.

## Structure
- `gene_panel_lists/`
  Versioned source inputs for each release. These folders contain gene lists and loci files.
- `adaptive_sampling_panels/`
  Versioned panel outputs for each release. Each release folder contains its own `README.md` describing what changed from the previous version.
- `scripts/`
  R scripts for the current three-step panel generation workflow.
- `run.sh`
  Example pipeline runner for a specific release configuration.
- `rproject.toml`, `rv.lock`, `rv/`
  Repo-local R environment configuration managed with `rv`.

## Using rv
This project uses `rv` to manage the R version, repositories, and package set locally for the repo.

Typical usage:

```sh
rv summary
rv sync
```

If you open an interactive R session from the repo root, the local [.Rprofile](/home/geonic/Documents/GitHub/gene-panels/.Rprofile) activates the repo `rv` environment automatically.

Useful commands:
- `rv summary`
  Show the current project configuration and status.
- `rv sync`
  Install or reconcile the project library with `rproject.toml` and `rv.lock`.
- `rv info --library`
  Print the repo-local library path used for this project.

The project repositories are defined in [rproject.toml](/home/geonic/Documents/GitHub/gene-panels/rproject.toml), and the repo-local package library under `rv/library/` is ignored by git.

## Current Workflow
The current panel generation workflow is:

1. `normalize_to_bed.R`
   Combines gene symbols and loci into a normalized minimal BED.
2. `quickpanel.R`
   Reads the normalized BED and applies padding to generate a sequencing panel BED/XLSX.
3. `merge_nearby_regions.R`
   Merges nearby padded intervals to avoid overlaps and writes the final merged BED/XLSX.

## Naming
By default, generated files follow:
- normalized: `YYYY-MM-DD-<version>-normalized-<padding_kb>kb.bed`
- panel: `YYYY-MM-DD-<version>-panel-<padding_kb>kb.bed/.xlsx`
- merged: `YYYY-MM-DD-<version>-merged.bed/.xlsx`

All three R scripts also accept `--output-file` to override the BED output path directly. For `quickpanel.R` and `merge_nearby_regions.R`, the Excel file path is derived from that BED path by swapping `.bed` to `.xlsx`.

## Release Tracking
- The panel output folders under `adaptive_sampling_panels/` mirror the release folders under `gene_panel_lists/`.
- Each release folder in `adaptive_sampling_panels/` documents the panel content differences relative to the previous version.
- Some older or manually named outputs may use legacy naming conventions; the per-release READMEs call those out where relevant.

## GitHub Releases
This repo includes a manual GitHub Actions workflow at [.github/workflows/publish-release-files.yml](/home/geonic/Documents/GitHub/gene-panels/.github/workflows/publish-release-files.yml) to publish selected files as release assets.

Use the `Publish Release Files` workflow from the Actions tab and provide:
- a `tag_name` such as `v1.0.5`
- a `release_name`
- newline-separated file paths or globs under `adaptive_sampling_panels/`

Example `files` input:

```text
adaptive_sampling_panels/v1.0.5-release/v1.0.5-analysis-set-0kb-pad.bed
adaptive_sampling_panels/v1.0.5-release/v1.0.5-analysis-set-0kb-pad.xlsx
adaptive_sampling_panels/v1.0.5-release/v1.0.5-sequencing-panel-20kb-pad-50kb-merge.bed
adaptive_sampling_panels/v1.0.5-release/v1.0.5-sequencing-panel-20kb-pad-50kb-merge.xlsx
```

The workflow creates the release if the tag does not exist yet, or updates the existing release and replaces matching assets if it does.

## Notes
To remove `0001_chr1_`-style prefixes from BED name fields for downstream tools such as MinKNOW coverage plotting:

```sh
awk 'BEGIN {OFS="\t"} {gsub(/[0-9]+_chr[0-9XY]+_/,"",$4); print}'
```
