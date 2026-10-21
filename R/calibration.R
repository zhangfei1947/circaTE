#' Build an empirical null calibrator for a declared simulation population
#'
#' Experimental API: never applies the development stress-mixture maps to
#' observed data. Users must attest all training datasets are entirely null,
#' specify the population context in te_data(), and use independent IDs for
#' evaluation. This checks metadata consistency, not the truth of that attestation.
#' Ties are inclusive, with add-one correction and no extrapolation below the
#' finite resolution floor. Scenario-max maps protect only the declared mixtures.
#'
#' @param null_fits A list of independent, uncalibrated simulation te_fit objects.
#' @param scenario One nonempty scenario label per fit.
#' @param all_null Must explicitly be TRUE to attest the training datasets are null.
#' @return A te_calibrator; no calibration is performed automatically by fit_te.
#' @export
build_te_calibrator <- function(null_fits, scenario, all_null = FALSE) {
  if (!isTRUE(all_null)) stop("Explicit all_null=TRUE attestation is required")
  if (!is.list(null_fits) || length(null_fits) < 2L ||
      !all(vapply(null_fits, inherits, logical(1), "te_fit"))) stop("Supply at least two null simulation fits")
  if (length(scenario) != length(null_fits) || anyNA(scenario) || any(!nzchar(as.character(scenario))))
    stop("Provide one nonempty scenario per null fit")
  first <- null_fits[[1]]; ids <- hashes <- character(length(null_fits)); groups <- list()
  for (i in seq_along(null_fits)) {
    f <- null_fits[[i]]
    .cte_calibration_ready(f)
    if (!identical(f$spec, first$spec) || !identical(f$versions, first$versions) ||
        !identical(f$data$calibration_context, first$data$calibration_context))
      stop("Calibration fits must share model specification, versions and simulation context")
    ids[i] <- f$data$dataset_id
    hashes[i] <- .cte_fingerprint(f$data)
    groups[[as.character(scenario[i])]] <- c(groups[[as.character(scenario[i])]], f$results$p)
  }
  if (anyDuplicated(ids)) stop("Calibration dataset IDs must be unique")
  if (anyDuplicated(hashes)) stop("Calibration inputs duplicate the same data under different IDs")
  groups <- lapply(groups, sort)
  structure(list(strata=groups, pooled=sort(unlist(groups,use.names=FALSE)),
    dataset_ids=ids, data_fingerprints=hashes, spec=first$spec, versions=first$versions,
    context=first$data$calibration_context, all_null_attested=TRUE,
    scope="Declared simulation mixture only; no real-data or subgroup conditional validity guarantee",
    created_utc=format(Sys.time(),tz="UTC",usetz=TRUE)), class="te_calibrator")
}
.cte_calibration_ready <- function(fit) {
  if (!inherits(fit, "te_fit") || fit$data$origin != "simulation")
    stop("Empirical calibration is restricted to simulation fits")
  if (!is.null(fit$calibration)) stop("Already calibrated fits cannot be used for calibration")
  validated <- te_data(fit$data$rna, fit$data$ribo, fit$data$metadata,
    offsets=list(rna=fit$data$off_rna,ribo=fit$data$off_ribo),origin=fit$data$origin,
    dataset_id=fit$data$dataset_id,calibration_context=fit$data$calibration_context)
  if (!identical(validated$meta,fit$data$meta)) stop("Fit metadata changed after fitting")
  if (is.null(fit$data$dataset_id) || is.null(fit$data$calibration_context))
    stop("Simulation calibration requires dataset_id and calibration_context")
  if (!is.data.frame(fit$results) || nrow(fit$results) != nrow(fit$data$rna) ||
      !identical(as.character(fit$results$gene_id),rownames(fit$data$rna)))
    stop("Calibration results must retain the complete aligned attempted gene family")
  if (!identical(fit$spec, .cte_spec(fit$data,fit$spec$method,fit$spec$period,fit$spec$trend)) ||
      !identical(fit$data_fingerprint,.cte_fingerprint(fit$data)))
    stop("Fit data or model specification changed after fitting")
  if (any(!is.finite(fit$results$p)) || any(fit$results$p < 0 | fit$results$p > 1))
    stop("Calibration requires finite p values in [0,1]")
}

#' Apply a frozen simulation-only null map
#' @param fit An independent, uncalibrated simulation te_fit.
#' @param calibrator A te_calibrator with matching specification and context.
#' @param mode Maximum cumulative null probability across scenarios, or pooled.
#' @return A te_fit with p_uncalibrated, calibrated p/q and explicit scope notes.
#' @export
calibrate_te <- function(fit, calibrator, mode = c("max", "pooled")) {
  mode <- match.arg(mode); .cte_calibration_ready(fit)
  if (!inherits(calibrator, "te_calibrator")) stop("calibrator must be te_calibrator")
  if (fit$data$dataset_id %in% calibrator$dataset_ids) stop("Evaluation dataset overlaps calibration training IDs")
  if (.cte_fingerprint(fit$data) %in% calibrator$data_fingerprints)
    stop("Evaluation counts overlap calibration training data despite different IDs")
  if (!identical(fit$spec, calibrator$spec) || !identical(fit$versions, calibrator$versions) ||
      !identical(fit$data$calibration_context, calibrator$context))
    stop("Fit does not match calibrator model, versions, design or simulation context")
  pp <- fit$results$p
  transform <- function(v) (1 + findInterval(pp, v))/(length(v)+1)
  p <- if (mode == "pooled") transform(calibrator$pooled) else Reduce(pmax,lapply(calibrator$strata,transform))
  fit$results$p_uncalibrated <- pp; fit$results$p <- p
  fit$results$q <- stats::p.adjust(p, "BH")
  floor <- if (mode=="pooled") 1/(length(calibrator$pooled)+1) else max(1/(lengths(calibrator$strata)+1))
  fit$calibration <- list(mode=mode, training_ids=calibrator$dataset_ids,
    context=calibrator$context, p_floor=floor, scope=calibrator$scope)
  fit$inferential_status <- "simulation_calibrated_experimental"
  fit$notes <- c(fit$notes, calibrator$scope, paste("Empirical p floor:",format(floor)))
  fit
}
