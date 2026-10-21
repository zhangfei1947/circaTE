# Internal paired count engines. Numerical designs and fitting options follow
# the frozen v1/v3 implementations; no simulation-derived calibration is applied.

.cte_count_require <- function(method) {
  if (method == "edger") {
    if (!requireNamespace("edgeR", quietly = TRUE) ||
        utils::packageVersion("edgeR") < "4.4.0")
      stop("The edger method requires edgeR >= 4.4.0.", call. = FALSE)
  } else if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Install the optional Bioconductor package DESeq2 to use this method.",
         call. = FALSE)
  }
}

.cte_count_design <- function(meta, period, trend, engine) {
  n <- nrow(meta)
  assay <- rep(c(0, 1), each = n)
  time <- rep(meta$time_hours, 2)
  phase <- 2 * pi * time / period
  if (engine == "edger") {
    pair <- rbind(diag(n), diag(n))
    colnames(pair) <- paste0("pair_", seq_len(n))
    full <- cbind(pair, assay_Ribo = assay)
    if (trend) full <- cbind(full, TE_trend = assay * (time - mean(time)) / period)
    full <- cbind(full, TE_cos = assay * cos(phase), TE_sin = assay * sin(phase))
    rr <- droplevels(factor(meta$replicate))
    if (nlevels(rr) > 1L) {
      r <- stats::model.matrix(~ rr)[, -1, drop = FALSE]
      r <- rbind(matrix(0, nrow(r), ncol(r)), r)
      colnames(r) <- paste0("TE_rep_", seq_len(ncol(r)))
      full <- cbind(full, r)
    }
    reduced <- full[, !colnames(full) %in% c("TE_cos", "TE_sin"), drop = FALSE]
  } else {
    # Preserve the original DESeq2 parameterization/order as well as its column space.
    full <- cbind(Intercept = 1, assay_ribo = assay)
    sample_factor <- factor(rep(seq_len(n), 2))
    mm <- stats::model.matrix(~ sample_factor)[, -1, drop = FALSE]
    colnames(mm) <- paste0("sample_", seq_len(ncol(mm)))
    full <- cbind(full, mm)
    if (trend) full <- cbind(full, ribo_trend = (time - mean(time)) / period * assay)
    rr <- droplevels(factor(rep(meta$replicate, 2)))
    if (nlevels(rr) > 1L) {
      z <- stats::model.matrix(~ rr)[, -1, drop = FALSE]
      colnames(z) <- paste0("ribo_rep_", seq_len(ncol(z)))
      full <- cbind(full, z * assay)
    }
    reduced <- full
    full <- cbind(full, ribo_cos = cos(phase) * assay, ribo_sin = sin(phase) * assay)
  }
  if (qr(full)$rank != ncol(full) || nrow(full) <= ncol(full))
    stop("The paired count design is rank deficient or has no residual degrees of freedom.",
         call. = FALSE)
  list(full = full, reduced = reduced, residual_df = nrow(full) - ncol(full))
}

.cte_count_edger <- function(counts, offsets, design, period) {
  y <- edgeR::DGEList(counts = counts)
  y$offset <- offsets
  warnings_seen <- character()
  fit <- withCallingHandlers({
    y <- edgeR::estimateDisp(y, design = design$full, robust = TRUE)
    edgeR::glmQLFit(y, design = design$full, robust = TRUE, legacy = FALSE)
  }, warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
  })
  test <- edgeR::glmQLFTest(fit,
    coef = match(c("TE_cos", "TE_sin"), colnames(design$full)))
  beta <- fit$coefficients[, c("TE_cos", "TE_sin"), drop = FALSE] / log(2)
  amp <- sqrt(rowSums(beta^2))
  phase <- (atan2(beta[, "TE_sin"], beta[, "TE_cos"]) %% (2 * pi)) * period / (2 * pi)
  phase[amp < sqrt(.Machine$double.eps)] <- NA_real_
  failed <- if (is.null(fit$failed)) rep(FALSE, nrow(counts)) else fit$failed
  failed <- failed | !is.finite(test$table$PValue) |
    !is.finite(beta[, 1]) | !is.finite(beta[, 2])
  out <- data.frame(gene_id = rownames(counts), p = test$table$PValue,
    cos = beta[, 1], sin = beta[, 2], amplitude_log2 = amp, phase_hours = phase,
    fit_failed = failed, informative = TRUE, row.names = NULL)
  attr(out, "diagnostics") <- list(engine = "edgeR", legacy = FALSE,
    package_version = as.character(utils::packageVersion("edgeR")),
    design_columns = colnames(design$full), design_rank = qr(design$full)$rank,
    residual_df = design$residual_df, test_df = 2L,
    NB_dispersion_scope = paste("Common NB baseline from highly expressed genes'",
      "dispersion trend; bias-adjusted genewise QL moderation"),
    NB_dispersion_summary = summary(fit$dispersion),
    QL_posterior_variance_summary = summary(fit$s2.post),
    adjusted_residual_df_summary = summary(fit$df.residual.adj),
    warnings = unique(warnings_seen),
    coefficient_estimation = "edgeR prior.count=0.125 shrinkage; QL test uses unshrunk fit")
  out
}

.cte_count_deseq <- function(counts, off_rna, off_ribo, design, period, finite) {
  sf <- exp(c(off_rna, off_ribo) - mean(c(off_rna, off_ribo)))
  if (any(!is.finite(sf)) || any(sf <= 0)) stop("Invalid DESeq2 size factors.", call. = FALSE)
  storage.mode(counts) <- "integer"
  n <- length(off_rna)
  colnames(counts) <- paste0(rep(c("RNA_", "RIBO_"), each = n), rep(seq_len(n), 2))
  cd <- data.frame(row.names = colnames(counts),
                   assay = factor(rep(c("RNA", "RIBO"), each = n)))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = cd,
                                      design = design$full)
  DESeq2::sizeFactors(dds) <- sf
  warnings_seen <- messages_seen <- character()
  dds <- withCallingHandlers(
    DESeq2::DESeq(dds, test = "LRT", full = design$full, reduced = design$reduced,
      betaPrior = FALSE, minReplicatesForReplace = Inf, parallel = FALSE, quiet = TRUE),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      messages_seen <<- c(messages_seen, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  res <- DESeq2::results(dds, cooksCutoff = FALSE, independentFiltering = FALSE)
  betas <- stats::coef(dds)
  mc <- as.data.frame(S4Vectors::mcols(dds))
  get_mc <- function(nm, default) if (nm %in% names(mc)) mc[[nm]] else rep(default, nrow(counts))
  full_ok <- get_mc("fullBetaConv", TRUE) & get_mc("betaConv", TRUE)
  reduced_ok <- get_mc("reducedBetaConv", TRUE)
  failed <- is.na(full_ok) | !full_ok | is.na(reduced_ok) | !reduced_ok |
    !is.finite(res$pvalue)
  bc <- betas[, "ribo_cos"]; bs <- betas[, "ribo_sin"]
  failed <- failed | !is.finite(bc) | !is.finite(bs)
  p <- res$pvalue
  if (finite) p <- stats::pf(pmax(res$stat, 0) / 2, df1 = 2,
                             df2 = design$residual_df, lower.tail = FALSE)
  failed <- failed | !is.finite(p)
  out <- data.frame(gene_id = rownames(counts), p = p, cos = bc, sin = bs,
    amplitude_log2 = sqrt(bc^2 + bs^2),
    phase_hours = (atan2(bs, bc) %% (2 * pi)) * period / (2 * pi),
    fit_failed = failed, informative = TRUE, lrt_statistic = res$stat,
    dispersion = DESeq2::dispersions(dds), converged_full = full_ok,
    converged_reduced = reduced_ok, p_raw = res$pvalue, row.names = NULL)
  attr(out, "diagnostics") <- list(engine = "DESeq2",
    package_version = as.character(utils::packageVersion("DESeq2")),
    design_columns = colnames(design$full), reduced_columns = colnames(design$reduced),
    design_rank = qr(design$full)$rank, residual_df = design$residual_df, test_df = 2L,
    calibration = if (finite) "heuristic_F_actual_residual_df" else "asymptotic_LRT",
    calibration_exact = FALSE, warnings = unique(warnings_seen),
    messages = unique(messages_seen),
    dispersion_fit_type = attr(DESeq2::dispersionFunction(dds), "fitType"),
    cooks_filter = FALSE, independent_filter = FALSE, outlier_replacement = FALSE,
    size_factors = sf)
  out
}

.cte_fit_count <- function(data, method = c("edger", "deseq2", "deseq2_finite"),
                           period = 24, trend = TRUE) {
  method <- match.arg(method)
  .cte_count_require(method)
  rna <- data$rna; ribo <- data$ribo; meta <- data$meta
  # Public validation is authoritative; these guards also protect direct internal use.
  if (!is.matrix(rna) || !is.matrix(ribo) || !identical(dim(rna), dim(ribo)) ||
      is.null(rownames(rna)) || !identical(rownames(rna), rownames(ribo)) ||
      ncol(rna) != nrow(meta)) stop("Aligned named count matrices and metadata are required.", call. = FALSE)
  counts <- cbind(rna, ribo)
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts)) ||
      any(counts > .Machine$integer.max)) stop("Finite nonnegative integer counts are required.", call. = FALSE)
  if (length(period) != 1L || !is.finite(period) || period <= 0)
    stop("period must be positive.", call. = FALSE)
  off <- c(data$off_rna, data$off_ribo)
  if (length(data$off_rna) != ncol(rna) || length(data$off_ribo) != ncol(rna) ||
      any(!is.finite(off))) stop("Finite per-sample log offsets are required.", call. = FALSE)
  informative <- rowSums(rna) > 0 & rowSums(ribo) > 0
  out <- data.frame(gene_id = rownames(rna), p = 1, cos = NA_real_, sin = NA_real_,
    amplitude_log2 = NA_real_, phase_hours = NA_real_, fit_failed = FALSE,
    informative = informative, row.names = NULL)
  d <- list(engine = if (method == "edger") "edgeR" else "DESeq2",
    no_informative_genes = !any(informative))
  if (any(informative)) {
    if (sum(informative) < 2L) stop("At least two informative genes are needed for count dispersion estimation.", call. = FALSE)
    engine <- if (method == "edger") "edger" else "deseq2"
    design <- .cte_count_design(meta, period, trend, engine)
    y <- counts[informative, , drop = FALSE]
    colnames(y) <- c(paste0("RNA_", seq_len(ncol(rna))), paste0("Ribo_", seq_len(ncol(ribo))))
    fit <- if (method == "edger") {
      offsets <- matrix(off, nrow(y), length(off), byrow = TRUE)
      .cte_count_edger(y, offsets, design, period)
    } else .cte_count_deseq(y, data$off_rna, data$off_ribo, design, period,
                            finite = method == "deseq2_finite")
    d <- attr(fit, "diagnostics")
    extra <- setdiff(names(fit), names(out))
    for (nm in extra) out[[nm]] <- NA
    out[informative, names(fit)] <- fit
    out$p[out$fit_failed] <- 1
    out[out$fit_failed, c("cos", "sin", "amplitude_log2", "phase_hours")] <- NA_real_
  }
  out$status <- ifelse(!out$informative, "uninformative",
                       ifelse(out$fit_failed, "fit_failed", "ok"))
  attr(out, "diagnostics") <- c(d, list(method = method, period_hours = period,
    trend = trend, paired = TRUE, serial_covariance_modeled = FALSE,
    genes_attempted = nrow(out), genes_informative = sum(out$informative),
    genes_uninformative = sum(!out$informative), genes_failed = sum(out$fit_failed),
    failed_gene_ids = out$gene_id[out$fit_failed],
    uninformative_policy = paste("Zero total in either assay: p=1, no effect estimate,",
      "excluded from dispersion learning, retained in testing family"),
    failure_policy = "Informative gene fitting failure: p=1 and fit_failed=TRUE"))
  out
}
