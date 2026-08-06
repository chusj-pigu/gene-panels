source("scripts/coordinate_sources.R")

expect_error <- function(expression, pattern) {
  error <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = function(condition) condition
  )

  stopifnot(inherits(error, "error"))
  stopifnot(grepl(pattern, conditionMessage(error), fixed = TRUE))
}

ucsc_called <- FALSE
ucsc_coordinates <- resolve_gene_coordinates(
  genes = "TP53",
  coordinate_source = "ucsc",
  genome = "hg38",
  fetchers = list(
    ucsc = function(genes, genome) {
      stopifnot(identical(genes, "TP53"), identical(genome, "hg38"))
      ucsc_called <<- TRUE
      data.frame(chrom = "chr17", start = 7661778L, end = 7687550L, label = "TP53")
    }
  )
)
stopifnot(ucsc_called)
stopifnot(identical(ucsc_coordinates$start, 7661778L))
stopifnot(identical(ucsc_coordinates$end, 7687550L))

parsed_ucsc_position <- parse_ucsc_browser_position(
  "chr17:7661779-7687546",
  label = "TP53"
)
stopifnot(identical(parsed_ucsc_position$chrom, "chr17"))
stopifnot(identical(parsed_ucsc_position$start, 7661778L))
stopifnot(identical(parsed_ucsc_position$end, 7687546L))

ensembl_coordinates <- resolve_gene_coordinates(
  genes = "TP53",
  coordinate_source = "ensembl",
  genome = "hg38",
  fetchers = list(
    ensembl = function(genes, genome) {
      data.frame(chrom = "17", start = 7661779L, end = 7687550L, label = "TP53")
    }
  )
)
stopifnot(identical(ensembl_coordinates$start, 7661778L))
stopifnot(identical(ensembl_coordinates$end, 7687550L))

expect_error(
  validate_coordinate_source("unknown"),
  "--coordinate-source must be one of"
)

expect_error(
  fetch_ensembl_gene_coordinates("TP53", genome = "hg19"),
  "Ensembl BioMart currently supports only --genome hg38"
)

expect_error(
  parse_ucsc_browser_position("not-a-genomic-position", label = "TP53"),
  "UCSC returned an unsupported genomic position"
)

cat("coordinate source tests passed\n")
