# Coordinate lookup helpers for normalize_to_bed.R.
#
# All returned coordinates use the 0-based, half-open convention required for
# BED output. Ensembl BioMart returns 1-based inclusive coordinates, whereas
# UCSC Genome Browser search coordinates are 1-based inclusive and are
# converted before they enter the normalized BED output.

coordinate_source_choices <- c("ensembl", "ucsc")

empty_gene_coordinates <- function() {
  data.frame(
    chrom = character(),
    start = integer(),
    end = integer(),
    label = character(),
    stringsAsFactors = FALSE
  )
}

validate_coordinate_source <- function(coordinate_source) {
  if (
    length(coordinate_source) != 1 ||
      is.na(coordinate_source) ||
      !coordinate_source %in% coordinate_source_choices
  ) {
    stop(
      "--coordinate-source must be one of: ",
      paste(coordinate_source_choices, collapse = ", "),
      call. = FALSE
    )
  }

  coordinate_source
}

coordinate_source_label <- function(coordinate_source) {
  switch(
    validate_coordinate_source(coordinate_source),
    ensembl = "Ensembl BioMart",
    ucsc = "UCSC Genome Browser REST API (HGNC track)"
  )
}

fetch_ensembl_gene_coordinates <- function(genes, genome) {
  if (!identical(genome, "hg38")) {
    stop(
      "Ensembl BioMart currently supports only --genome hg38 in this workflow; ",
      "use --coordinate-source ucsc for another UCSC assembly.",
      call. = FALSE
    )
  }

  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  records <- biomaRt::getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
    filters = "hgnc_symbol",
    values = genes,
    mart = ensembl
  )

  data.frame(
    chrom = records$chromosome_name,
    start = records$start_position,
    end = records$end_position,
    label = records$hgnc_symbol,
    stringsAsFactors = FALSE
  )
}

fetch_ucsc_gene_coordinates <- function(genes, genome) {
  if (!requireNamespace("curl", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "The UCSC coordinate source requires curl and jsonlite. Run `rv sync` from the repository root.",
      call. = FALSE
    )
  }

  records <- lapply(seq_along(genes), function(index) {
    if (index < length(genes)) {
      on.exit(Sys.sleep(0.1), add = TRUE)
    }

    gene <- genes[[index]]
    query_url <- paste0(
      "https://api.genome.ucsc.edu/search?search=",
      utils::URLencode(gene, reserved = TRUE),
      "&genome=",
      utils::URLencode(genome, reserved = TRUE)
    )
    response <- curl::curl_fetch_memory(query_url)

    if (response$status_code != 200L) {
      stop(
        "UCSC API request failed for ", gene,
        " on ", genome,
        " with HTTP status ", response$status_code,
        call. = FALSE
      )
    }

    payload <- jsonlite::fromJSON(rawToChar(response$content), simplifyVector = FALSE)
    position_matches <- payload$positionMatches
    if (is.null(position_matches)) {
      return(empty_gene_coordinates())
    }

    hgnc_tracks <- Filter(
      function(track) identical(track$name, "hgnc"),
      position_matches
    )
    hgnc_matches <- unlist(
      lapply(hgnc_tracks, function(track) track$matches),
      recursive = FALSE
    )
    exact_matches <- Filter(
      function(match) {
        identical(toupper(as.character(match$posName)), toupper(gene))
      },
      hgnc_matches
    )

    if (length(exact_matches) == 0) {
      return(empty_gene_coordinates())
    }

    intervals <- lapply(exact_matches, function(match) {
      parse_ucsc_browser_position(match$position, label = gene)
    })
    do.call(rbind, intervals)
  })

  records <- records[vapply(records, nrow, integer(1)) > 0]
  if (length(records) == 0) {
    return(empty_gene_coordinates())
  }

  do.call(rbind, records)
}

parse_ucsc_browser_position <- function(position, label) {
  fields <- regmatches(
    position,
    regexec("^([^:]+):([0-9]+)-([0-9]+)$", position)
  )[[1]]

  if (length(fields) != 4) {
    stop("UCSC returned an unsupported genomic position: ", position, call. = FALSE)
  }

  start <- as.integer(fields[[3]])
  end <- as.integer(fields[[4]])
  if (is.na(start) || is.na(end) || start < 1L || end < start) {
    stop("UCSC returned an invalid genomic position: ", position, call. = FALSE)
  }

  data.frame(
    chrom = fields[[2]],
    start = start - 1L,
    end = end,
    label = label,
    stringsAsFactors = FALSE
  )
}

normalize_gene_coordinates <- function(records, coordinate_source, genes) {
  coordinate_source <- validate_coordinate_source(coordinate_source)
  required_columns <- c("chrom", "start", "end", "label")
  missing_columns <- setdiff(required_columns, names(records))
  if (length(missing_columns) > 0) {
    stop(
      "Coordinate source returned incomplete data; missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  coordinates <- data.frame(
    chrom = as.character(records$chrom),
    start = suppressWarnings(as.integer(records$start)),
    end = suppressWarnings(as.integer(records$end)),
    label = as.character(records$label),
    stringsAsFactors = FALSE
  )

  coordinates <- coordinates[
    !is.na(coordinates$label) & coordinates$label %in% genes,
    ,
    drop = FALSE
  ]

  if (nrow(coordinates) == 0) {
    return(empty_gene_coordinates())
  }

  if (any(is.na(coordinates$start)) || any(is.na(coordinates$end))) {
    stop("Coordinate source returned non-numeric start or end positions.", call. = FALSE)
  }

  if (identical(coordinate_source, "ensembl")) {
    if (any(coordinates$start < 1L)) {
      stop("Ensembl returned a start position below 1.", call. = FALSE)
    }
    coordinates$start <- coordinates$start - 1L
  }

  if (any(coordinates$start < 0L) || any(coordinates$end <= coordinates$start)) {
    stop("Coordinate source returned invalid genomic intervals.", call. = FALSE)
  }

  coordinates
}

resolve_gene_coordinates <- function(
  genes,
  coordinate_source,
  genome,
  fetchers = list(
    ensembl = fetch_ensembl_gene_coordinates,
    ucsc = fetch_ucsc_gene_coordinates
  )
) {
  coordinate_source <- validate_coordinate_source(coordinate_source)
  if (length(genes) == 0) {
    return(empty_gene_coordinates())
  }

  fetcher <- fetchers[[coordinate_source]]
  if (!is.function(fetcher)) {
    stop("No fetcher is configured for coordinate source: ", coordinate_source, call. = FALSE)
  }

  normalize_gene_coordinates(
    fetcher(genes, genome),
    coordinate_source = coordinate_source,
    genes = genes
  )
}
