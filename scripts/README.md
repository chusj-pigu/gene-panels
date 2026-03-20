# Scripts

The `scripts/` folder contains the current three-step adaptive sampling panel workflow.

## Pipeline
1. `normalize_to_bed.R`
2. `quickpanel.R`
3. `merge_nearby_regions.R`

## normalize_to_bed.R
Builds a minimal BED from a gene list, a loci file, or both.

### Inputs
- `--genes`
  Text file with one gene symbol per line.
- `--loci`
  BED-like tab-delimited file with at least `chrom`, `start`, `end`, and optional `name`.

### Behavior
- queries Ensembl BioMart for gene coordinates
- warns about duplicated input gene symbols and removes them before querying
- warns about gene symbols with no BioMart hit
- warns about post-normalization duplicate names and removes them by default
- can restrict output to canonical chromosomes with `--canonical-only`

### Output
- default BED name: `YYYY-MM-DD-<version>-normalized-<padding_kb>kb.bed`
- override with `--output-file`

## quickpanel.R
Takes a normalized BED-like file and applies padding to produce a sequencing panel BED and Excel workbook.

### Inputs
- `--input`
  Normalized BED-like file from `normalize_to_bed.R`
- `--padding`
  Padding in base pairs

### Output
- default BED/XLSX names: `YYYY-MM-DD-<version>-panel-<padding_kb>kb.*`
- override the BED path with `--output-file`
- the XLSX path is derived from the BED path

## merge_nearby_regions.R
Merges nearby padded intervals to collapse overlapping or near-overlapping targets.

### Inputs
- `--input`
  BED-like panel file, typically from `quickpanel.R`
- `--gap`
  Maximum allowed gap in bp between intervals to merge

### Behavior
- prints a merge report listing which source regions were collapsed
- merged region IDs include the contributing gene/region names
- writes an Excel workbook with:
  - `Panel`
  - `Merged_Genes`
  - `Summary`

### Output
- default BED/XLSX names: `YYYY-MM-DD-<version>-merged.*`
- override the BED path with `--output-file`
- the XLSX path is derived from the BED path

## Shared Conventions
- all scripts require `--version`
- all scripts create `--outdir` if needed
- all scripts accept `--output-file`; if omitted, they fall back to their default filename constructor
- `quickpanel.R` and `merge_nearby_regions.R` derive the Excel path from the BED output path

## Example
```sh
Rscript scripts/normalize_to_bed.R \
  --version v1.0.5 \
  --genes gene_panel_lists/v1.0.5-20260316-release/v1.0.5-20260316-release-genes.txt \
  --loci gene_panel_lists/v1.0.5-20260316-release/v1.0.5-20260316-release-loci.txt \
  --padding 0 \
  --outdir adaptive_sampling_panels/v1.0.5-20260316-release

Rscript scripts/quickpanel.R \
  --input adaptive_sampling_panels/v1.0.5-20260316-release/YYYY-MM-DD-v1.0.5-normalized-0kb.bed \
  --version v1.0.5 \
  --padding 20000 \
  --outdir adaptive_sampling_panels/v1.0.5-20260316-release

Rscript scripts/merge_nearby_regions.R \
  --input adaptive_sampling_panels/v1.0.5-20260316-release/YYYY-MM-DD-v1.0.5-panel-20kb.bed \
  --version v1.0.5 \
  --gap 50000 \
  --outdir adaptive_sampling_panels/v1.0.5-20260316-release
```
