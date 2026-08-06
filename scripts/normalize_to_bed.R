#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
})

script_file <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
)
source(file.path(dirname(normalizePath(script_file)), "coordinate_sources.R"))

# ---- Define options ----
option_list <- list(
  make_option(c("-g", "--genes"), type = "character", default = NULL,
              help = "Path to the gene list file", metavar = "FILE"),
  make_option(c("-l", "--loci"), type = "character", default = NULL,
              help = "Optional loci file", metavar = "FILE"),
  make_option(c("-v", "--version"), type = "character", default = NULL,
              help = "Version token used in output naming", metavar = "STRING"),
  make_option(c("-o", "--output-file"), type = "character", default = NULL, dest = "output_file",
              help = "Explicit output BED file path; overrides generated naming", metavar = "FILE"),
  make_option(c("-p", "--padding"), type = "integer", default = 0,
              help = "Padding for start and end positions", metavar = "NUMBER"),
  make_option(c("--coordinate-source"), type = "character", default = "ensembl", dest = "coordinate_source",
              help = "Gene coordinate source: ensembl (BioMart) or ucsc (UCSC HGNC track)", metavar = "SOURCE"),
  make_option(c("--genome"), type = "character", default = "hg38",
              help = "Genome assembly; Ensembl currently supports hg38, UCSC uses this assembly directly", metavar = "ASSEMBLY"),
  make_option(c("-r", "--outdir"), type = "character", default = ".",
              help = "Output directory", metavar = "DIR"),
  make_option(c("-H", "--header"), action = "store_true", default = FALSE,
              help = "Loci file includes a header row"),
  make_option(c("-c", "--canonical-only"), action = "store_true", default = FALSE, dest = "canonical_only",
              help = "Restrict output to chromosomes 1-22, X, and Y"),
  make_option(c("--keep-post-normalization-duplicates"), action = "store_false", default = TRUE,
              dest = "remove_post_normalization_duplicates",
              help = "Keep post-normalization duplicate entries instead of removing them")
)

# ---- Parse options ----
opt <- parse_args(OptionParser(option_list = option_list))
cat("📁 Output directory:", opt$outdir, "\n")
cat("🏷️ Version:", opt$version, "\n")
cat("📄 Output file override:", ifelse(is.null(opt$output_file), "<generated>", opt$output_file), "\n")
cat("🧬 Padding:", opt$padding, "\n")
cat("🧭 Coordinate source:", opt$coordinate_source, "\n")
cat("🧬 Genome assembly:", opt$genome, "\n")
cat("📚 Header row present:", isTRUE(opt$header), "\n")
cat("🧪 Canonical only:", isTRUE(opt$canonical_only), "\n")
cat("🧹 Remove post-normalization duplicates:", isTRUE(opt$remove_post_normalization_duplicates), "\n")

# ---- Validate options ----
if (is.null(opt$genes) && is.null(opt$loci)) {
  stop("❌ Provide at least one of --genes or --loci.")
}

if (is.null(opt$version) || !nzchar(opt$version)) {
  stop("❌ Provide a non-empty --version.")
}

if (opt$padding < 0) {
  stop("❌ --padding must be zero or greater.")
}

opt$coordinate_source <- validate_coordinate_source(opt$coordinate_source)

if (is.null(opt$genome) || !nzchar(opt$genome)) {
  stop("❌ Provide a non-empty --genome.")
}

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
}

# ---- Prepare gene list ----
genes <- character(0)
if (!is.null(opt$genes)) {
  if (!file.exists(opt$genes)) {
    stop("❌ Gene list file not found: ", opt$genes)
  }

  cat("📖 Reading gene list from:", opt$genes, "\n")
  genes <- readLines(opt$genes) |>
    str_trim() |>
    discard(~ .x == "")

  duplicated_genes <- genes[duplicated(genes)] |> unique()
  if (length(duplicated_genes) > 0) {
    warning(
      paste0(
        "Removing ",
        length(duplicated_genes),
        " duplicated gene symbol(s): ",
        paste(duplicated_genes, collapse = ", ")
      ),
      call. = FALSE
    )
    genes <- unique(genes)
  }

  cat("🧬 Loaded", length(genes), "genes\n")
}

# ---- Helpers ----
normalize_chromosome <- function(x) {
  x |>
    as.character() |>
    str_trim() |>
    str_remove("^chr") |>
    str_to_upper() |>
    str_replace("^M$", "MT")
}

format_ucsc_chromosome <- function(x) {
  if_else(x == "MT", "chrM", paste0("chr", x))
}

build_region_id <- function(label, chrom, start, end, fallback_prefix = "region") {
  cleaned_label <- label |>
    replace_na("") |>
    str_trim()

  fallback <- paste0(fallback_prefix, "_chr", chrom, "_", start, "_", end)
  if_else(cleaned_label == "", fallback, cleaned_label)
}

canonical_levels <- c(as.character(1:22), "X", "Y")

# ---- Retrieve gene loci from the selected coordinate source ----
gene_regions <- tibble()
if (length(genes) > 0) {
  source_label <- coordinate_source_label(opt$coordinate_source)
  cat("🔌 Connecting to", source_label, "...\n")
  cat("🔍 Querying", source_label, "for gene coordinates...\n")
  gene_regions <- resolve_gene_coordinates(
    genes = genes,
    coordinate_source = opt$coordinate_source,
    genome = opt$genome
  ) |>
    as_tibble() |>
    transmute(
      source = "gene",
      chrom = normalize_chromosome(chrom),
      start = as.integer(start),
      end = as.integer(end),
      label = as.character(label)
    )

  cat("✅ Retrieved", nrow(gene_regions), "gene loci from", source_label, "\n")

  matched_symbols <- unique(gene_regions$label)
  missing_symbols <- setdiff(genes, matched_symbols)

  if (length(missing_symbols) > 0) {
    warning(
      paste0(
        "No ", source_label, " hit for ",
        length(missing_symbols),
        " gene symbol(s): ",
        paste(missing_symbols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

# ---- Load loci file ----
loci_regions <- tibble()
if (!is.null(opt$loci)) {
  if (!file.exists(opt$loci)) {
    stop("❌ Loci file not found: ", opt$loci)
  }

  cat("➕ Reading loci from:", opt$loci, "\n")
  loci_raw <- read_tsv(
    opt$loci,
    col_names = isTRUE(opt$header),
    comment = "#",
    show_col_types = FALSE,
    progress = FALSE
  )

  if (ncol(loci_raw) < 3) {
    stop("❌ Loci file must have at least 3 columns: chrom, start, end.")
  }

  loci_label <- if (ncol(loci_raw) >= 4) {
    as.character(loci_raw[[4]])
  } else {
    rep(NA_character_, nrow(loci_raw))
  }

  loci_regions <- loci_raw |>
    as_tibble() |>
    rename(
      chrom = 1,
      start = 2,
      end = 3
    ) |>
    mutate(
      source = "loci",
      label = loci_label,
      chrom = normalize_chromosome(chrom),
      start = suppressWarnings(as.integer(start)),
      end = suppressWarnings(as.integer(end))
    ) |>
    select(source, chrom, start, end, label)

  cat("📄 Loaded", nrow(loci_regions), "loci rows\n")
}

# ---- Combine and normalize ----
regions <- bind_rows(gene_regions, loci_regions)
initial_region_count <- nrow(regions)

if (nrow(regions) == 0) {
  stop("❌ No regions were collected from the provided inputs.")
}

if (any(is.na(regions$start)) || any(is.na(regions$end))) {
  stop("❌ Start and end positions must be numeric in all inputs.")
}

if (any(regions$end < regions$start)) {
  stop("❌ Found regions where end < start.")
}

non_canonical_regions <- tibble()
if (isTRUE(opt$canonical_only)) {
  non_canonical_regions <- regions |>
    filter(!chrom %in% canonical_levels)
  regions <- regions |>
    filter(chrom %in% canonical_levels)

  cat("🧪 Retained", nrow(regions), "canonical regions\n")
  if (nrow(non_canonical_regions) > 0) {
    cat("🗑️ Dropped", nrow(non_canonical_regions), "non-canonical region(s):\n")
    print(non_canonical_regions)
  }
}

if (nrow(regions) == 0) {
  stop("❌ No regions remain after filtering.")
}

chrom_levels <- c(
  canonical_levels,
  "MT",
  sort(setdiff(unique(regions$chrom), c(canonical_levels, "MT")))
)

cat("📦 Normalizing regions to minimal BED...\n")
bed_regions_raw <- regions |>
  mutate(
    start = pmax(0L, start - opt$padding),
    end = end + opt$padding,
    chrom = factor(chrom, levels = chrom_levels)
  ) |>
  arrange(chrom, start, end, label) |>
  mutate(
    chrom = as.character(chrom),
    name = build_region_id(label, chrom, start, end),
    chrom = format_ucsc_chromosome(chrom)
  ) |>
  select(chrom, start, end, name)

exact_duplicate_entries <- bed_regions_raw |>
  add_count(chrom, start, end, name, name = "entry_count") |>
  filter(entry_count > 1) |>
  arrange(name, chrom, start, end)

if (nrow(exact_duplicate_entries) > 0) {
  cat(
    "🗑️ Removed ",
    nrow(exact_duplicate_entries),
    " exact duplicate normalized entr",
    ifelse(nrow(exact_duplicate_entries) == 1, "y", "ies"),
    ":\n",
    sep = ""
  )
  print(exact_duplicate_entries)
}

bed_regions <- bed_regions_raw |>
  distinct(chrom, start, end, name, .keep_all = TRUE)

mitochondrial_region_count <- sum(bed_regions$chrom == "chrM")
if (mitochondrial_region_count > 0) {
  cat("🧬 Rewrote", mitochondrial_region_count, "mitochondrial region(s) to UCSC chromosome name chrM\n")
}

duplicate_normalized_names <- bed_regions |>
  add_count(name, name = "name_count") |>
  filter(name_count > 1) |>
  arrange(name, chrom, start, end)

if (nrow(duplicate_normalized_names) > 0) {
  duplicated_symbols <- unique(duplicate_normalized_names$name)
  warning(
    paste0(
      "Post-normalization duplicate gene symbol(s) remain for ",
      length(duplicated_symbols),
      " symbol(s): ",
      paste(duplicated_symbols, collapse = ", ")
    ),
    call. = FALSE
  )
  cat("📝 Post-normalization duplicate entries:\n")
  print(duplicate_normalized_names)

  if (isTRUE(opt$remove_post_normalization_duplicates)) {
    cat("⚠️ Removing post-normalization duplicates by symbol can drop valid multi-locus targets if the same gene name appears at distinct loci\n")
    duplicate_entries_removed <- duplicate_normalized_names |>
      group_by(name) |>
      slice(-1) |>
      ungroup()
    cat(
      "🗑️ Removed ",
      nrow(duplicate_entries_removed),
      " post-normalization duplicate entr",
      ifelse(nrow(duplicate_entries_removed) == 1, "y", "ies"),
      ":\n",
      sep = ""
    )
    print(duplicate_entries_removed)
    bed_regions <- bed_regions |>
      arrange(name, chrom, start, end) |>
      distinct(name, .keep_all = TRUE)

    cat("🧹 Removed post-normalization duplicate entries, keeping the first entry per symbol\n")
  }
}

cat("✅ Normalized", nrow(bed_regions), "regions\n")
cat("📉 Normalization summary:\n")
cat("  - Started with", initial_region_count, "region row(s)\n")
if (isTRUE(opt$canonical_only)) {
  cat("  - Dropped", nrow(non_canonical_regions), "non-canonical region row(s)\n")
}
cat("  - Removed", nrow(exact_duplicate_entries), "exact duplicate normalized row(s)\n")
if (nrow(duplicate_normalized_names) > 0 && isTRUE(opt$remove_post_normalization_duplicates)) {
  duplicate_rows_removed_count <- duplicate_normalized_names |>
    group_by(name) |>
    slice(-1) |>
    ungroup() |>
    nrow()
  cat("  - Removed", duplicate_rows_removed_count, "duplicate-symbol row(s) after normalization\n")
}
cat("  - Kept", nrow(bed_regions), "final region row(s)\n")

summary_df <- tibble(
  Input_Genes = length(genes),
  Gene_Regions = nrow(gene_regions),
  Loci_Regions = nrow(loci_regions),
  Output_Regions = nrow(bed_regions),
  Padding_bp = opt$padding
)

# ---- Write output BED ----
basename <- paste0(Sys.Date(), "-", opt$version, "-normalized-", opt$padding / 1000, "kb")
output_path <- if (is.null(opt$output_file)) {
  file.path(opt$outdir, paste0(basename, ".bed"))
} else {
  opt$output_file
}

cat("💾 Writing BED...\n")
bed_regions |>
  write_tsv(output_path, col_names = FALSE)

cat("✅ Minimal BED written to:", output_path, "\n")
print(summary_df)
