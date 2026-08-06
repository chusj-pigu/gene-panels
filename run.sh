#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TODAY="$(date +%F)"

GENOME="hg38"
COORDINATE_SOURCE="${COORDINATE_SOURCE:-ensembl}"
PADDING_BP=20000
MERGE_GAP_BP=50000
OUTDIR="${ROOT_DIR}/adaptive_sampling_panels/v1.0.5-release"

GENESPATH="${ROOT_DIR}/gene_panel_lists/v1.0.5-20260316-release/v1.0.5-20260316-release-genes.txt"
LOCIPATH="${ROOT_DIR}/gene_panel_lists/v1.0.5-20260316-release/v1.0.5-20260316-release-loci.txt"

VERSION="v1.0.5"

NORMALIZED_BED="${OUTDIR}/${VERSION}-normalized-inputs.bed"
ANALYSIS_BED="${OUTDIR}/${VERSION}-analysis-set-0kb-pad.bed"
PREMERGE_BED="${OUTDIR}/${VERSION}-pre-merge-panel-$((PADDING_BP / 1000))kb-pad.bed"
SEQUENCING_BED="${OUTDIR}/${VERSION}-sequencing-panel-$((PADDING_BP / 1000))kb-pad-$((MERGE_GAP_BP / 1000))kb-merge.bed"

Rscript "${ROOT_DIR}/scripts/normalize_to_bed.R" \
    --canonical-only \
    --version="${VERSION}" \
    --genes="${GENESPATH}" \
    --loci="${LOCIPATH}" \
    --coordinate-source="${COORDINATE_SOURCE}" \
    --genome="${GENOME}" \
    --padding=0 \
    --outdir="${OUTDIR}" \
    --output-file="${NORMALIZED_BED}"

Rscript "${ROOT_DIR}/scripts/quickpanel.R" \
    --input="${NORMALIZED_BED}" \
    --version="${VERSION}" \
    --padding=0 \
    --genome="${GENOME}" \
    --outdir="${OUTDIR}" \
    --output-file="${ANALYSIS_BED}"

Rscript "${ROOT_DIR}/scripts/quickpanel.R" \
    --input="${NORMALIZED_BED}" \
    --version="${VERSION}" \
    --padding="${PADDING_BP}" \
    --genome="${GENOME}" \
    --outdir="${OUTDIR}" \
    --output-file="${PREMERGE_BED}"

Rscript "${ROOT_DIR}/scripts/merge_nearby_regions.R" \
    --input="${PREMERGE_BED}" \
    --version="${VERSION}" \
    --gap="${MERGE_GAP_BP}" \
    --genome="${GENOME}" \
    --outdir="${OUTDIR}" \
    --output-file="${SEQUENCING_BED}"
