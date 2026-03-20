#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
})

# ---- Define options ----
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to the normalized BED-like regions file", metavar = "FILE"),
  make_option(c("-v", "--version"), type = "character", default = NULL,
              help = "Version token used in output naming", metavar = "STRING"),
  make_option(c("-o", "--output-file"), type = "character", default = NULL, dest = "output_file",
              help = "Explicit output BED file path; overrides generated naming", metavar = "FILE"),
  make_option(c("-p", "--padding"), type = "integer", default = 20000,
              help = "Padding for start and end positions", metavar = "NUMBER"),
  make_option(c("-d", "--genome"), type = "character", default = "hg38",
              help = "UCSC genome version", metavar = "STRING"),
  make_option(c("-H", "--header"), action = "store_true", default = FALSE,
              help = "Input file includes a header row"),
  make_option(c("-r", "--outdir"), type = "character", default = ".",
              help = "Output directory", metavar = "DIR")
)

# ---- Parse options ----
opt <- parse_args(OptionParser(option_list = option_list))
cat("📁 Output directory:", opt$outdir, "\n")
cat("🏷️ Version:", opt$version, "\n")
cat("📄 Output file override:", ifelse(is.null(opt$output_file), "<generated>", opt$output_file), "\n")
cat("🧬 Padding:", opt$padding, "\n")
cat("🧬 Genome:", opt$genome, "\n")
cat("📚 Header row present:", isTRUE(opt$header), "\n")

# ---- Validate inputs ----
if (is.null(opt$input) || !file.exists(opt$input)) {
  stop("❌ Please provide a valid normalized BED-like file with --input.")
}

if (is.null(opt$version) || !nzchar(opt$version)) {
  stop("❌ Provide a non-empty --version.")
}

if (opt$padding < 0) {
  stop("❌ --padding must be zero or greater.")
}

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
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

build_region_id <- function(label, chrom, row_number, start, end) {
  cleaned_label <- label |>
    replace_na("") |>
    str_trim()

  region_index <- str_pad(row_number, 4, pad = "0")
  fallback <- paste0("region_", region_index, "_chr", chrom, "_", start, "_", end)
  final_label <- if_else(cleaned_label == "", fallback, cleaned_label)

  paste0(region_index, "_chr", chrom, "_", final_label)
}

canonical_levels <- c(as.character(1:22), "X", "Y")

# ---- Read normalized regions ----
cat("📖 Reading normalized regions from:", opt$input, "\n")
regions_raw <- read_tsv(
  opt$input,
  col_names = isTRUE(opt$header),
  comment = "#",
  show_col_types = FALSE,
  progress = FALSE
)

if (ncol(regions_raw) < 3) {
  stop("❌ Input must contain at least 3 tab-delimited columns: chrom, start, end.")
}

region_name <- if (ncol(regions_raw) >= 4) {
  as.character(regions_raw[[4]])
} else {
  rep(NA_character_, nrow(regions_raw))
}

regions <- regions_raw |>
  as_tibble() |>
  rename(
    chrom = 1,
    start = 2,
    end = 3
  ) |>
  mutate(
    region_name = region_name,
    chrom = normalize_chromosome(chrom),
    start = suppressWarnings(as.integer(start)),
    end = suppressWarnings(as.integer(end))
  ) |>
  select(chrom, start, end, region_name)

if (any(is.na(regions$start)) || any(is.na(regions$end))) {
  stop("❌ Start and end positions must be numeric.")
}

if (any(regions$end < regions$start)) {
  stop("❌ Found regions where end < start.")
}

chrom_levels <- c(
  canonical_levels,
  "MT",
  sort(setdiff(unique(regions$chrom), c(canonical_levels, "MT")))
)

cat("📦 Generating padded panel data...\n")
panel_regions <- regions |>
  mutate(chrom = factor(chrom, levels = chrom_levels)) |>
  arrange(chrom, start, end, region_name) |>
  mutate(
    chrom = as.character(chrom),
    start = pmax(0L, start - opt$padding),
    end = end + opt$padding,
    row_number = row_number(),
    gene_id = build_region_id(region_name, chrom, row_number, start, end),
    strand = ".",
    chrom = format_ucsc_chromosome(chrom),
    size = end - start
  ) |>
  select(chrom, start, end, gene_id, strand, size) |>
  mutate(
    ucsc_link = paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=", opt$genome,
                       "&position=", chrom, ":", start, "-", end),
    ensembl_link = paste0("https://www.ensembl.org/Homo_sapiens/Location/View?r=",
                          chrom, ":", start, "-", end)
  )

mitochondrial_region_count <- sum(panel_regions$chrom == "chrM")
if (mitochondrial_region_count > 0) {
  cat("🧬 Rewrote", mitochondrial_region_count, "mitochondrial region(s) to UCSC chromosome name chrM\n")
}

cat("✅ Prepared", nrow(panel_regions), "panel regions\n")

# ---- Coverage stats ----
total_genome_size <- ifelse(opt$genome == "hg38", 3.2e9, 3.09e9)
panel_size <- sum(panel_regions$end - panel_regions$start)
percent_coverage <- (panel_size / total_genome_size) * 100
cat("📏 Panel covers", round(percent_coverage, 4), "% of genome\n")

summary_df <- tibble(
  Input_Regions = nrow(regions),
  Panel_Regions = nrow(panel_regions),
  Padding_bp = opt$padding,
  Total_Genome_Size = total_genome_size,
  Panel_Size = panel_size,
  Percent_Genome_Coverage = percent_coverage
)

# ---- Write output files ----
basename <- paste0(Sys.Date(), "-", opt$version, "-panel-", opt$padding / 1000, "kb")
bed_output_path <- if (is.null(opt$output_file)) {
  file.path(opt$outdir, paste0(basename, ".bed"))
} else {
  opt$output_file
}
xlsx_output_path <- sub("\\.bed$", ".xlsx", bed_output_path)
if (identical(xlsx_output_path, bed_output_path)) {
  xlsx_output_path <- paste0(bed_output_path, ".xlsx")
}

cat("💾 Writing Excel summary and BED...\n")
writexl::write_xlsx(
  list(Panel = panel_regions, Summary = summary_df),
  xlsx_output_path
)

panel_regions |>
  select(chrom, start, end, gene_id) |>
  write_tsv(bed_output_path, col_names = FALSE)

cat("✅ BED and Excel written to:", opt$outdir, "\n")
print(summary_df)
