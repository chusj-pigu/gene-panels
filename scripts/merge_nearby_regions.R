#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
})

build_merged_region_id <- function(chrom, row_number, source_regions) {
  region_index <- str_pad(row_number, 4, pad = "0")
  cleaned_sources <- source_regions |>
    replace_na("") |>
    str_split(";") |>
    map_chr(function(parts) {
      if (length(parts) == 0) {
        return("")
      }

      parts |>
        str_trim() |>
        str_replace("^\\d+_chr[^_]+_", "") |>
        paste(collapse = "_")
    }) |>
    str_replace_all(";", "_") |>
    str_replace_all("[^[:alnum:]_]+", "_") |>
    str_replace_all("_+", "_") |>
    str_replace("^_", "") |>
    str_replace("_$", "")

  suffix <- if_else(cleaned_sources == "", "merged", cleaned_sources)
  paste0(region_index, "_chr", chrom, "_", suffix)
}

# ---- Define options ----
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Path to the BED-like regions file", metavar = "FILE"),
  make_option(c("-v", "--version"), type = "character", default = NULL,
              help = "Version token used in output naming", metavar = "STRING"),
  make_option(c("-o", "--output-file"), type = "character", default = NULL, dest = "output_file",
              help = "Explicit output BED file path; overrides generated naming", metavar = "FILE"),
  make_option(c("-g", "--gap"), type = "integer", default = 0,
              help = "Maximum gap in bp allowed between regions to merge", metavar = "NUMBER"),
  make_option(c("-d", "--genome"), type = "character", default = "hg38",
              help = "UCSC genome version", metavar = "STRING"),
  make_option(c("-r", "--outdir"), type = "character", default = ".",
              help = "Output directory", metavar = "DIR"),
  make_option(c("-H", "--header"), action = "store_true", default = FALSE,
              help = "Input file includes a header row")
)

# ---- Parse options ----
opt <- parse_args(OptionParser(option_list = option_list))
cat("📁 Output directory:", opt$outdir, "\n")
cat("🏷️ Version:", opt$version, "\n")
cat("📄 Output file override:", ifelse(is.null(opt$output_file), "<generated>", opt$output_file), "\n")
cat("🧬 Merge gap:", opt$gap, "bp\n")
cat("🧬 Genome:", opt$genome, "\n")
cat("📚 Header row present:", isTRUE(opt$header), "\n")

# ---- Validate inputs ----
if (is.null(opt$input) || !file.exists(opt$input)) {
  stop("❌ Please provide a valid regions file with --input.")
}

if (is.null(opt$version) || !nzchar(opt$version)) {
  stop("❌ Provide a non-empty --version.")
}

if (opt$gap < 0) {
  stop("❌ --gap must be zero or greater.")
}

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
}

# ---- Read regions ----
cat("📖 Reading regions from:", opt$input, "\n")
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

regions <- regions_raw |>
  as_tibble() |>
  rename(
    chrom = 1,
    start = 2,
    end = 3
  ) |>
  mutate(
    region_label = if (ncol(regions_raw) >= 4) as.character(.data[[names(regions_raw)[4]]]) else NA_character_,
    chrom = as.character(chrom),
    start = suppressWarnings(as.integer(start)),
    end = suppressWarnings(as.integer(end))
  )

if (any(is.na(regions$start)) || any(is.na(regions$end))) {
  stop("❌ Start and end columns must be numeric.")
}

if (any(regions$end < regions$start)) {
  stop("❌ Found regions where end < start.")
}

cat("🧬 Loaded", nrow(regions), "regions\n")

# ---- Normalize and sort ----
chrom_levels <- c(as.character(1:22), "X", "Y", "M", "MT")
chrom_levels <- c(
  chrom_levels,
  sort(setdiff(unique(str_remove(regions$chrom, "^chr")), chrom_levels))
)

regions_sorted <- regions |>
  mutate(chrom = str_remove(chrom, "^chr")) |>
  mutate(chrom = factor(chrom, levels = chrom_levels)) |>
  arrange(chrom, start, end) |>
  mutate(chrom = as.character(chrom))

cat("🧪 Retained", nrow(regions_sorted), "regions after normalization\n")

# ---- Merge nearby regions ----
cat("🔗 Merging nearby regions...\n")
merged_regions <- regions_sorted |>
  group_by(chrom) |>
  mutate(
    running_end = lag(cummax(end)),
    starts_new_group = if_else(is.na(running_end) | start > running_end + opt$gap, 1L, 0L),
    merge_group = cumsum(starts_new_group)
  ) |>
  group_by(chrom, merge_group) |>
  summarise(
    start = min(start),
    end = max(end),
    input_region_count = n(),
    merged_size = max(end) - min(start),
    source_regions = paste(
      unique(na.omit(region_label)),
      collapse = ";"
    ),
    .groups = "drop"
  ) |>
  arrange(factor(chrom, levels = chrom_levels), start, end) |>
  mutate(
    row_number = row_number(),
    region_id = build_merged_region_id(chrom, row_number, source_regions)
  ) |>
  mutate(
    chrom = paste0("chr", chrom),
    size = end - start,
    source_regions = na_if(source_regions, "")
  ) |>
  mutate(
    ucsc_link = paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=", opt$genome,
                       "&position=", chrom, ":", start, "-", end),
    ensembl_link = paste0("https://www.ensembl.org/Homo_sapiens/Location/View?r=",
                          chrom, ":", start, "-", end)
  ) |>
  select(chrom, start, end, region_id, input_region_count, merged_size, size, source_regions, ucsc_link, ensembl_link)

cat("✅ Reduced", nrow(regions_sorted), "regions to", nrow(merged_regions), "merged regions\n")

merged_genes_df <- merged_regions |>
  transmute(
    region_id,
    chrom,
    start,
    end,
    source_regions = replace_na(source_regions, "")
  ) |>
  separate_rows(source_regions, sep = ";", convert = FALSE) |>
  mutate(source_region = na_if(source_regions, "")) |>
  select(region_id, chrom, start, end, source_region)

# ---- Print merge report ----
merge_report <- merged_regions |>
  filter(input_region_count > 1) |>
  transmute(
    report_line = paste0(
      region_id,
      " -> ",
      chrom,
      ":",
      start,
      "-",
      end,
      " | merged ",
      input_region_count,
      " regions",
      if_else(
        is.na(source_regions),
        "",
        paste0(" | sources: ", source_regions)
      )
    )
  )

if (nrow(merge_report) > 0) {
  cat("📝 Merged groups:\n")
  cat(paste0("  - ", merge_report$report_line, collapse = "\n"), "\n")
} else {
  cat("📝 Merged groups: none\n")
}

# ---- Coverage stats ----
original_span <- sum(regions_sorted$end - regions_sorted$start)
merged_span <- sum(merged_regions$end - merged_regions$start)
net_span_change <- merged_span - original_span
total_genome_size <- ifelse(opt$genome == "hg38", 3.2e9, 3.09e9)
percent_coverage <- (merged_span / total_genome_size) * 100

summary_df <- tibble(
  Input_Regions = nrow(regions_sorted),
  Merged_Regions = nrow(merged_regions),
  Region_Count_Reduction = nrow(regions_sorted) - nrow(merged_regions),
  Merge_Gap_bp = opt$gap,
  Total_Genome_Size = total_genome_size,
  Input_Total_bp = original_span,
  Merged_Total_bp = merged_span,
  Net_Span_Change_bp = net_span_change,
  Percent_Genome_Coverage = percent_coverage
)

cat("📏 Input span:", original_span, "bp\n")
cat("📏 Merged span:", merged_span, "bp\n")
cat("📐 Net span change:", net_span_change, "bp\n")
cat("📏 Panel covers", round(percent_coverage, 4), "% of genome\n")

# ---- Write output files ----
basename <- paste0(Sys.Date(), "-", opt$version, "-merged")
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
  list(Panel = merged_regions, Merged_Genes = merged_genes_df, Summary = summary_df),
  xlsx_output_path
)

merged_regions |>
  select(chrom, start, end, region_id) |>
  write_tsv(bed_output_path, col_names = FALSE)

cat("✅ BED and Excel written to:", opt$outdir, "\n")
print(summary_df)
