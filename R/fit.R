#' Fit a fixed-period paired translational-efficiency model
#'
#' The default joint edgeR model absorbs shared RNA changes with sample-pair
#' fixed effects. Its reference distribution does not model serial dependence.
#' GLS uses estimated observation precision and culture-block AR1 plus nugget;
#' plug-in covariance uncertainty and low-count log bias remain. All methods
#' return exploratory conditional inference, not universally calibrated FDR.
#' Single-assay-zero genes receive p=1 and stay in BH; no low-count observations
#' are removed. Fits never select the smallest p across methods or periods.
#'
#' @param data A te_data object, one biological group.
#' @param method edger, weighted_ar (moderated), unweighted_ar (matched moderated),
#'   weighted_independent (moderated), weighted_ar_ordinary, deseq2, or deseq2_finite.
#' @param period One pre-specified positive period in hours; default 24.
#' @param trend Include an assay-specific linear trend nuisance term.
#' @return A te_fit containing results, diagnostics, input, and model specification.
#'   amplitude_log2 is the harmonic half-range; phase_hours is the harmonic peak
#'   modulo period relative to time zero, excluding the linear trend.
#' @export
#' @examples
#' d <- example_te_data()
#' f <- fit_te(d)
#' head(f$results)
fit_te <- function(data, method = c("edger", "weighted_ar", "unweighted_ar",
                    "weighted_independent", "weighted_ar_ordinary", "deseq2", "deseq2_finite"),
                   period = 24, trend = TRUE) {
  if (!inherits(data, "te_data")) stop("data must be created with te_data()")
  method <- match.arg(method)
  if (!is.numeric(period) || length(period) != 1L || !is.finite(period) || period <= 0)
    stop("period must be one positive finite number")
  if (!is.logical(trend) || length(trend) != 1L || is.na(trend)) stop("trend must be TRUE or FALSE")
  # Revalidate public objects after possible manual modification, preserving offsets.
  data <- te_data(data$rna, data$ribo, data$metadata,
    offsets = list(rna = data$off_rna, ribo = data$off_ribo), origin = data$origin,
    dataset_id = data$dataset_id, calibration_context = data$calibration_context) |>
    .cte_restore_normalization(data)
  n <- ncol(data$rna); b <- factor(data$meta$replicate)
  X <- matrix(1, n, 1)
  if (trend) X <- cbind(X, (data$meta$time_hours - mean(data$meta$time_hours)) / period)
  if (nlevels(b) > 1) X <- cbind(X, stats::model.matrix(~ b)[, -1, drop = FALSE])
  ph <- 2*pi*data$meta$time_hours/period
  X <- cbind(X, cos(ph), sin(ph))
  singular <- svd(X, nu=0L, nv=0L)$d
  if (qr(X)$rank != ncol(X) || sum(singular > max(singular)*1e-8) != ncol(X) || n <= ncol(X))
    stop("TE design is rank deficient, numerically aliased, or has no residual degrees of freedom")
  ar <- method %in% c("weighted_ar", "unweighted_ar", "weighted_ar_ordinary")
  if (ar && nlevels(b) < 2L) stop("AR models require at least two independent culture batches")
  notes <- c("Exploratory conditional model p-values; universal FDR control is not established.",
             "Relative TE from sequencing; neither absolute translation nor proof of endogenous circadian rhythm.")
  if (diff(range(data$meta$time_hours)) < 2*period)
    notes <- c(notes, "Observed span covers fewer than two specified periods.")
  if (method %in% c("edger", "deseq2", "deseq2_finite"))
    tab <- .cte_fit_count(data, method, period, trend)
  else tab <- .cte_fit_gls(data, method, period, trend)
  diagnostics <- attr(tab, "diagnostics")
  if (!identical(as.character(tab$gene_id), rownames(data$rna))) stop("Internal result gene alignment failure")
  bad <- !is.finite(tab$p) | tab$p < 0 | tab$p > 1
  if (!"fit_failed" %in% names(tab)) tab$fit_failed <- bad
  tab$fit_failed <- tab$fit_failed | bad
  informative <- rowSums(data$rna) > 0 & rowSums(data$ribo) > 0
  tab$informative <- informative
  tab$p[bad | !informative] <- 1
  tab$fit_failed[!informative] <- FALSE
  tab[!informative,c("cos","sin","amplitude_log2","phase_hours")] <- NA_real_
  tab$q <- stats::p.adjust(tab$p, "BH")
  tab$method <- method
  tab$status <- ifelse(!informative, "uninformative", ifelse(tab$fit_failed, "fit_failed", "fitted"))
  if (!is.null(diagnostics$converged) && !isTRUE(diagnostics$converged)) {
    warning("Covariance optimization did not converge; inspect diagnostics", call. = FALSE)
    notes <- c(notes, "Covariance optimization did not converge.")
  }
  if (isTRUE(diagnostics$boundary)) notes <- c(notes, "Covariance estimate is near an optimization boundary.")
  if (any(tab$fit_failed)) notes <- c(notes, "Failed gene tests retained at p=1 in the attempted family.")
  spec <- .cte_spec(data, method, period, trend)
  structure(list(results = tab, diagnostics = diagnostics, data = data, spec = spec,
    notes = unique(notes), inferential_status = "exploratory_conditional",
    versions = .cte_versions(), data_fingerprint = .cte_fingerprint(data)), class = "te_fit")
}
.cte_restore_normalization <- function(validated, original) {
  validated$normalization <- original$normalization
  validated$normalization_factors <- original$normalization_factors
  validated
}
.cte_versions <- function() {
  pkgs <- c("circaTE", "edgeR", "limma", "statmod", "DESeq2")
  setNames(vapply(pkgs, function(p) if (requireNamespace(p, quietly=TRUE))
    as.character(utils::packageVersion(p)) else NA_character_, character(1)), pkgs)
}

#' @export
print.te_fit <- function(x, ...) {
  cat("circaTE", x$spec$method, "|", x$spec$period, "hours |", nrow(x$results), "genes\n")
  cat(sum(x$results$q < .05), "nominal BH q<0.05; status:", x$inferential_status, "\n")
  cat(paste0("- ", x$notes, collapse="\n"), "\n")
  invisible(x)
}

.cte_spec <- function(data, method, period, trend) {
  list(method = method, period = period, trend = trend, harmonics = 1L,
    normalization = data$normalization, n_genes = nrow(data$rna), n_pairs = ncol(data$rna),
    time_hours = data$meta$time_hours,
    batch_structure = as.integer(factor(data$meta$replicate)),
    package_version = as.character(utils::packageVersion("circaTE")))
}
.cte_fingerprint <- function(data) {
  # Exclude dataset labels so a renamed copy cannot enter the evaluation split.
  i <- order(rownames(data$rna)); j <- order(data$metadata$sample_id)
  digest::digest(list(rna=unname(data$rna[i,j,drop=FALSE]),
    ribo=unname(data$ribo[i,j,drop=FALSE]), time=data$meta$time_hours[j],
    batch=as.integer(factor(data$meta$replicate[j])),
    off_rna=unname(data$off_rna[j]),off_ribo=unname(data$off_ribo[j])),algo="sha256")
}
