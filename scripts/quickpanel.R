#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    library(optparse)
    library(biomaRt)
    library(tidyverse)
    library(dplyr)
})
# ---- Define options ----
option_list <- list(
    make_option(c("-g", "--genes"), type = "character", default = "sj-panel/genes_list.txt",
                            help = "Path to the gene list file [default: %default]", metavar = "FILE"),
    make_option(c("-o", "--output"), type = "character", default = "sj-panel/CHUSJ-2025-Adaptive-Panel-LP-Genes.bed",
                            help = "Path to the output file [default: %default]", metavar = "FILE"),
    make_option(c("-p", "--padding"), type = "integer", default = 20000,
                            help = "Padding for start and end positions [default: %default]", metavar = "NUMBER"),
    make_option(c("-d", "--genome"), type = "character", default = "hg38",
                            help = "UCSC genome version (e.g., hg38, hg19) [default: %default]", metavar = "STRING"),
    make_option(c("-l", "--loci"), type = "character", default = NULL,
                            help = "Path to the file with additional loci [default: %default]", metavar = "FILE"),
    make_option(c("-r", "--outdir"), type = "character", default = ".",
                            help = "Output directory [default: %default]", metavar = "DIR")
)

# ---- Parse options ----
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Read the gene list ----
genes <- readLines(opt$genes)

# ---- Connect to Ensembl BioMart ----
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# ---- Retrieve genomic loci information ----
gene_loci <- getBM(
    attributes = c(
        "hgnc_symbol",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand"
    ),
    filters = "hgnc_symbol",
    values = genes,
    mart = ensembl
)

# ---- Filter and format the data ----
gene_loci_canonical <- gene_loci |>
    filter(chromosome_name %in% c(as.character(1:22), "X", "Y"))

# ---- Add manual loci ----
# Parse only if loci not empty
if( !is.null(opt$loci) ) {
    read_tsv(opt$loci, col_types = c("cciii")) |>
        bind_rows(gene_loci_canonical) |>
        mutate(chromosome_name = factor(chromosome_name, levels = c(as.character(1:22), "X", "Y"))) ->
        gene_loci_canonical
}

# ---- Add padding, generate gene IDs, and prepare for BED format ----
gene_loci_bed <- gene_loci_canonical |>
    arrange(chromosome_name, start_position) |>
    mutate(row_number = row_number()) |>
    mutate(gene_index = str_pad(row_number, 4, pad = "0")) |>
    mutate(gene_id = paste0(gene_index, "_", "chr", chromosome_name, "_", hgnc_symbol)) |>
    mutate(strand = if_else(strand == 1, "+", "-")) |>
    rename(start = start_position, end = end_position, chrom = chromosome_name) |>
    mutate(start = pmax(0, start - opt$padding), end = end + opt$padding) |>
    mutate(chrom = paste0("chr", chrom)) |>
    mutate(size = end - start) |>
    select(chrom, start, end, gene_id, strand, size) |>
    mutate(
        ucsc_link = paste0(
            "https://genome.ucsc.edu/cgi-bin/hgTracks?db=",
            opt$genome,
            "&position=",
            chrom,
            ":",
            start,
            "-",
            end
        )
    ) |>
    mutate(
        ensembl_link = paste0(
            "https://www.ensembl.org/Homo_sapiens/Location/View?r=",
            chrom,
            ":",
            start,
            "-",
            end
        )
    )

# ---- Calculate genome coverage ----
total_genome_size <- ifelse(opt$genome == "hg38", 3.2e9, 3090000000) # Approximate size of human genome
panel_size <- sum(gene_loci_bed$end - gene_loci_bed$start)
percent_coverage <- (panel_size / total_genome_size) * 100

# ---- Create a summary dataframe ----
summary_df <- tibble(
    Total_Genome_Size = total_genome_size,
    Panel_Size = panel_size,
    Percent_Genome_Coverage = percent_coverage
)

# ---- Add summary dataframe to xlsx list ----

gene_loci_out <- 
    list(
        "Panel" = gene_loci_bed, 
        "Summary" = summary_df
        )

gene_loci_out |>
    writexl::write_xlsx(paste0(opt$outdir,"/",Sys.Date(), "-CHUSJ-Panel-",opt$padding/1000, "kb-",  opt$output, ".xlsx"))


gene_loci_bed |>
    select(chrom, start, end, gene_id) |>
    write_tsv(paste0(opt$outdir,"/", Sys.Date(), "-CHUSJ-Panel-",opt$padding/1000, "kb-", opt$output, ".bed"), col_names = FALSE)

# ---- Print summary ----
print(summary_df)
