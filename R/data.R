#' Validate and prepare paired RNA and Ribo counts
#'
#' One object represents one biological group. Both matrices must use identical
#' gene rows and biological sample columns, in metadata order. No implicit
#' reordering, expression filtering, or imputation is performed. Supply full
#' count matrices for normalization, or externally computed log offsets if the
#' matrices were already subsetted. Offsets are natural log effective library
#' sizes, not raw sizes. Zero-count genes remain in the attempted BH family.
#'
#' @param rna,ribo Numeric integer count matrices, genes by paired samples.
#' @param metadata Data frame with sample_id, time_hours, batch; optional group
#'   must have one distinct value. batch identifies a culture followed over time.
#' @param normalization Either separate-assay TMM or library_size.
#' @param offsets Optional list with rna and ribo numeric log-offset vectors.
#' @param origin Observed data or a simulation. Calibration accepts simulation only.
#' @param dataset_id Unique dataset identifier, required for simulation calibration.
#' @param calibration_context User-defined simulation population identifier.
#' @return A te_data object with aligned counts, metadata, offsets and provenance.
#' @export
te_data <- function(rna, ribo, metadata, normalization = c("TMM", "library_size"),
                    offsets = NULL, origin = c("observed", "simulation"),
                    dataset_id = NULL, calibration_context = NULL) {
  normalization <- match.arg(normalization); origin <- match.arg(origin)
  check_counts <- function(x, assay) {
    if (!is.matrix(x) || !is.numeric(x) || length(dim(x)) != 2L ||
        nrow(x) < 2L || ncol(x) < 2L)
      stop(assay, " must be a numeric matrix with >=2 genes and >=2 samples")
    if (any(!is.finite(x)) || any(x < 0) || any(x != round(x)))
      stop(assay, " counts must be finite nonnegative integers")
    if (any(x > .Machine$integer.max)) stop(assay, " counts exceed supported integer range")
    for (ids in dimnames(x)) {
      if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
        stop(assay, " requires unique nonempty gene and sample names")
    }
    if (is.null(rownames(x)) || is.null(colnames(x))) stop(assay, " requires dimension names")
    if (any(colSums(x) <= 0)) stop(assay, " contains an empty library")
  }
  check_counts(rna, "RNA"); check_counts(ribo, "Ribo")
  if (!identical(dimnames(rna), dimnames(ribo)))
    stop("RNA and Ribo gene/sample names and order must match exactly")
  if (!is.data.frame(metadata) || nrow(metadata) != ncol(rna) ||
      !all(c("sample_id", "time_hours", "batch") %in% names(metadata)))
    stop("metadata must match samples and contain sample_id, time_hours, batch")
  ids <- as.character(metadata$sample_id)
  if (!identical(ids, colnames(rna))) stop("metadata sample_id order must match count columns exactly")
  if (!is.numeric(metadata$time_hours) || any(!is.finite(metadata$time_hours)))
    stop("time_hours must be finite numeric hours")
  b <- as.character(metadata$batch)
  if (anyNA(b) || any(!nzchar(b))) stop("batch labels must be complete and nonempty")
  if (anyDuplicated(data.frame(batch = b, time = metadata$time_hours)))
    stop("Each batch must have unique time points; do not label unrelated replicates as one culture")
  if ("group" %in% names(metadata) &&
      (anyNA(metadata$group) || length(unique(metadata$group)) != 1L))
    stop("One te_data object must contain exactly one biological group")
  if (length(unique(metadata$time_hours)) < 3L) stop("At least three distinct times are required")
  for (v in list(dataset_id, calibration_context))
    if (!is.null(v) && (!is.character(v) || length(v) != 1L || is.na(v) || !nzchar(v)))
      stop("dataset_id and calibration_context must be nonempty single strings")
  metadata$sample_id <- ids; metadata$batch <- b
  meta <- data.frame(bio_sample_id = ids, time_hours = metadata$time_hours, replicate = b)
  factors <- list()
  if (is.null(offsets)) {
    get_offset <- function(x) {
      d <- edgeR::DGEList(counts = x)
      if (normalization == "TMM") d <- edgeR::calcNormFactors(d, method = "TMM")
      list(offset = log(d$samples$lib.size * d$samples$norm.factors),
           factors = d$samples$norm.factors)
    }
    a <- get_offset(rna); z <- get_offset(ribo)
    off_rna <- a$offset; off_ribo <- z$offset
    factors <- list(rna = a$factors, ribo = z$factors)
    normalization_used <- normalization
  } else {
    if (!is.list(offsets) || !all(c("rna", "ribo") %in% names(offsets)))
      stop("offsets must be a list with rna and ribo log vectors")
    for (v in offsets[c("rna", "ribo")]) {
      if (!is.numeric(v) || !is.null(dim(v)) || length(v) != ncol(rna) ||
          any(!is.finite(v)) || any(!is.finite(exp(v))) || any(exp(v) <= 0))
        stop("offsets must be finite log effective-library-size vectors")
      if (!is.null(names(v)) && !identical(names(v), ids)) stop("Named offsets must follow sample_id order")
    }
    off_rna <- unname(offsets$rna); off_ribo <- unname(offsets$ribo)
    normalization_used <- "supplied_log_offsets"
  }
  structure(list(rna = rna, ribo = ribo, metadata = metadata, meta = meta,
    off_rna = off_rna, off_ribo = off_ribo, normalization = normalization_used,
    normalization_factors = factors, origin = origin, dataset_id = dataset_id,
    calibration_context = calibration_context), class = "te_data")
}

#' Read tab-separated paired count files
#' @param rna_file,ribo_file TSV files with gene_id first and sample IDs as headers.
#' @param metadata_file TSV file containing sample_id, time_hours and batch.
#' @param ... Passed to te_data.
#' @return A te_data object.
#' @export
read_te_data <- function(rna_file, ribo_file, metadata_file, ...) {
  read_counts <- function(path) {
    d <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!ncol(d) || names(d)[1] != "gene_id") stop("First count column must be gene_id")
    ids <- as.character(d[[1]]); d <- d[-1]
    if (!all(vapply(d, is.numeric, logical(1)))) stop("Count columns must be numeric")
    m <- as.matrix(d); rownames(m) <- ids; m
  }
  te_data(read_counts(rna_file), read_counts(ribo_file),
          utils::read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE), ...)
}

#' @export
print.te_data <- function(x, ...) {
  cat("circaTE paired data:", nrow(x$rna), "genes;", ncol(x$rna), "pairs;",
      length(unique(x$meta$replicate)), "culture batches\n")
  cat("Normalization:", x$normalization, "| origin:", x$origin, "\n")
  invisible(x)
}
