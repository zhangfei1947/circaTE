# Internal paired logTE GLS engine, ported from frozen v2/v3 implementations.
# Covariance and precision are estimated; reference p values are uncalibrated.
# Numerical lineage: circadian_te_lowcount_v3_20261021/methods/weighted_gls.R.

.cte_ar_prepare <- function(rna, ribo, meta, off_rna, off_ribo, period, trend, block_rep) {
  rna <- as.matrix(rna); ribo <- as.matrix(ribo)
  if (!is.numeric(rna) || !is.numeric(ribo) || !identical(dim(rna), dim(ribo)) ||
      is.null(rownames(rna)) || !identical(rownames(rna), rownames(ribo)) ||
      anyDuplicated(rownames(rna)) || any(!is.finite(c(rna, ribo))) ||
      any(c(rna, ribo) < 0)) stop("Require aligned nonnegative count matrices and unique gene names")
  n <- ncol(rna)
  if (nrow(meta) != n || !is.numeric(meta$time_hours) ||
      any(!is.finite(meta$time_hours)) || is.null(meta$replicate) || anyNA(meta$replicate))
    stop("Require matched metadata with finite numeric time_hours and culture replicate")
  if (length(period) != 1 || !is.finite(period) || period <= 0) stop("Invalid period")
  if (length(off_rna) != n || length(off_ribo) != n ||
      any(!is.finite(c(off_rna, off_ribo)))) stop("Require one finite offset per sample")
  lr <- exp(off_rna); lb <- exp(off_ribo)
  if (any(!is.finite(c(lr, lb))) || any(c(lr, lb) <= 0)) stop("Invalid effective library sizes")
  te <- log2(sweep(ribo + 0.5, 2, lb + 1, "/")) -
    log2(sweep(rna + 0.5, 2, lr + 1, "/"))
  b <- factor(meta$replicate)
  if (anyDuplicated(paste(b, meta$time_hours, sep = ":")))
    stop("Each culture batch must have unique time points")
  dt <- diff(sort(unique(meta$time_hours)))
  if (!length(dt)) stop("At least two distinct times are required")
  step <- min(dt)
  lag <- abs(outer(meta$time_hours, meta$time_hours, "-")) / step
  regular <- max(abs(lag - round(lag))) < 1e-7
  if (regular) lag <- round(lag)
  same <- outer(as.character(b), as.character(b), "==")
  X0 <- matrix(1, n, 1, dimnames = list(NULL, "intercept"))
  if (trend) X0 <- cbind(X0, time = (meta$time_hours - mean(meta$time_hours)) / period)
  if (block_rep && nlevels(b) > 1) X0 <- cbind(X0, stats::model.matrix(~ b)[, -1, drop = FALSE])
  w <- 2 * pi * meta$time_hours / period
  X1 <- cbind(X0, cos = cos(w), sin = sin(w))
  X2 <- cbind(X1, cos2 = cos(2 * w), sin2 = sin(2 * w))
  # Package engines estimate and test H1. Unused H2 columns must not reject an
  # identifiable H1 design (for example, a four-phase six-hour schedule).
  if (qr(X1)$rank < ncol(X1) || n <= ncol(X1)) stop("H1 design needs full rank and residual df")
  list(y = t(te), genes = rownames(rna), X0 = X0, X1 = X1, X2 = X2,
       lag = lag, same = same, step = step, regular = regular, period = period)
}

.cte_ar_covariance <- function(x, rho, nugget = 0) {
  if (length(rho) != 1 || !is.finite(rho) || abs(rho) >= 1 ||
      length(nugget) != 1 || !is.finite(nugget) || nugget < 0 || nugget > 1)
    stop("Require |rho|<1 and 0<=nugget<=1")
  if (rho < 0 && !x$regular) stop("Negative rho requires a regular time lattice")
  R <- (1 - nugget) * rho^x$lag * x$same
  diag(R) <- 1
  R
}


.cte_wg_prepare <- function(rna, ribo, meta, off_rna, off_ribo,
                        period = 24, trend = TRUE, block_rep = TRUE) {
  x <- .cte_ar_prepare(rna, ribo, meta, off_rna, off_ribo, period, trend, block_rep)
  if (!requireNamespace("limma", quietly = TRUE)) stop("limma required")
  # Saturated time means accommodate RNA troughs of arbitrary shape. The assay
  # mean model estimates precisions only; it is not the TE rhythm test design.
  tt <- factor(meta$time_hours); bb <- factor(meta$replicate)
  W <- if (block_rep && nlevels(bb) > 1) stats::model.matrix(~ tt + bb) else stats::model.matrix(~ tt)
  if (qr(W)$rank != ncol(W) || nrow(W) <= ncol(W)) stop("Precision design needs residual df")
  a <- limma::voom(rna, design = W, lib.size = exp(off_rna), normalize.method = "none", plot = FALSE)
  b <- limma::voom(ribo, design = W, lib.size = exp(off_ribo), normalize.method = "none", plot = FALSE)
  v <- 1 / a$weights + 1 / b$weights
  if (any(!is.finite(v)) || any(v <= 0)) stop("Invalid voom variances")
  x$precision <- t(1 / v)
  x$informative <- rowSums(rna) > 0 & rowSums(ribo) > 0
  x$precision_design <- colnames(W)
  x$precision_residual_df <- nrow(W) - ncol(W)
  x
}

# Batched Cholesky over genes: arrays contain one p x p matrix per gene.
# D differs for every gene. Multiplication by R^-1 is shared, but design
# whitening is NOT shared; it is recomputed from D_g^-1 X for every gene.
.cte_wg_gls <- function(x, X, rho, nugget) {
  R <- .cte_ar_covariance(x, rho, nugget)
  Ri <- chol2inv(chol(R)); sw <- sqrt(x$precision)
  yy <- x$y * sw; iy <- Ri %*% yy
  p <- ncol(X); ng <- ncol(yy)
  sx <- lapply(seq_len(p), function(j) sw * X[, j])
  ix <- lapply(sx, function(z) Ri %*% z)
  A <- array(0, c(ng, p, p)); rhs <- matrix(0, ng, p)
  for (j in seq_len(p)) {
    rhs[, j] <- colSums(sx[[j]] * iy)
    for (k in seq_len(j)) A[, j, k] <- A[, k, j] <- colSums(sx[[j]] * ix[[k]])
  }
  L <- array(0, dim(A))
  for (j in seq_len(p)) for (k in seq_len(j)) {
    z <- A[, j, k]
    if (k > 1) for (h in seq_len(k - 1)) z <- z - L[, j, h] * L[, k, h]
    if (j == k) {
      if (any(z <= 0) || any(!is.finite(z))) stop("Ill-conditioned weighted design")
      L[, j, k] <- sqrt(z)
    } else L[, j, k] <- z / L[, k, k]
  }
  u <- beta <- matrix(0, ng, p)
  for (j in seq_len(p)) {
    z <- rhs[, j]
    if (j > 1) for (k in seq_len(j - 1)) z <- z - L[, j, k] * u[, k]
    u[, j] <- z / L[, j, j]
  }
  for (j in rev(seq_len(p))) {
    z <- u[, j]
    if (j < p) for (k in (j + 1):p) z <- z - L[, k, j] * beta[, k]
    beta[, j] <- z / L[, j, j]
  }
  # Compute residuals explicitly instead of subtracting two large sums.
  residual <- (x$y - X %*% t(beta)) * sw
  rss <- pmax(colSums(residual * (Ri %*% residual)), .Machine$double.xmin)
  logdetX <- Reduce(`+`, lapply(seq_len(p), function(j) 2 * log(L[, j, j])))
  colnames(beta) <- colnames(X)
  list(beta = beta, rss = rss, logdetX = logdetX,
       logdetR = as.numeric(determinant(R, logarithm = TRUE)$modulus))
}

.cte_wg_estimate <- function(x) {
  df <- nrow(x$X1) - ncol(x$X1)
  usable <- colSums(stats::lm.fit(x$X1, x$y)$residuals^2) > 1e-14
  if (!is.null(x$informative)) usable <- usable & x$informative
  if (!any(usable)) return(list(rho = 0, nugget = 1, converged = TRUE,
    boundary = TRUE, objective = NA_real_, genes = 0))
  xx <- x; xx$y <- x$y[, usable, drop = FALSE]; xx$precision <- x$precision[, usable, drop = FALSE]
  bounds <- if (x$regular) c(-.6, .95) else c(0, .95)
  # log|D_g| is constant in rho/nugget, so omitted from the objective.
  fn <- function(par) {
    z <- .cte_wg_gls(xx, xx$X1, par[1], par[2])
    sum(z$logdetR + z$logdetX + df * log(z$rss / df))
  }
  opts <- lapply(list(c(.3, .2), c(.75, .6)), function(start)
    stats::optim(start, fn, method = "L-BFGS-B", lower = c(bounds[1], 0),
      upper = c(bounds[2], .98), control = list(factr = 1e8, maxit = 80)))
  best <- opts[[which.min(vapply(opts, function(z) z$value, numeric(1)))]]
  list(rho = best$par[1], nugget = best$par[2], converged = best$convergence == 0,
    boundary = min(abs(best$par[1] - bounds)) < .002 || best$par[2] < .002 || best$par[2] > .978,
    objective = best$value, genes = sum(usable), bounds = bounds)
}

.cte_wg_test <- function(x, rho, nugget, moderate = FALSE) {
  full <- .cte_wg_gls(x, x$X1, rho, nugget)
  null <- .cte_wg_gls(x, x$X0, rho, nugget)
  df <- nrow(x$X1) - ncol(x$X1); q <- ncol(x$X1) - ncol(x$X0)
  s2 <- full$rss / df; df2 <- rep(df, length(s2))
  if (moderate) {
    eb <- limma::squeezeVar(s2, df = df, robust = TRUE)
    s2 <- eb$var.post; df2 <- pmin(df + eb$df.prior, length(s2) * df)
  }
  Fstat <- pmax(null$rss - full$rss, 0) / q / s2
  p <- stats::pf(Fstat, q, df2, lower.tail = FALSE)
  constant <- apply(x$y, 2, function(y) diff(range(y)) < 1e-12)
  if (!is.null(x$informative)) constant <- constant | !x$informative
  p[constant] <- 1; Fstat[constant] <- 0
  co <- full$beta; co[constant, c("cos", "sin")] <- 0
  amp <- sqrt(co[, "cos"]^2 + co[, "sin"]^2)
  phase <- (atan2(co[, "sin"], co[, "cos"]) %% (2*pi)) * x$period / (2*pi)
  phase[amp < 1e-12] <- NA_real_
  data.frame(gene_id = x$genes, p = p, cos = co[, "cos"], sin = co[, "sin"],
    amplitude_log2 = amp, phase_hours = phase, F = Fstat, residual_df = df,
    denominator_df = df2, test_df = q, rho = rho, nugget = nugget,
    mean_precision = colMeans(x$precision), min_precision = apply(x$precision, 2, min),
    max_precision = apply(x$precision, 2, max),
    informative = if (is.null(x$informative)) TRUE else x$informative, row.names = NULL)
}

.cte_wg_annotate <- function(ans, x, est, method, moderate) {
  ans$method <- method
  attr(ans, "diagnostics") <- c(est, list(
    precision_design = x$precision_design, precision_residual_df = x$precision_residual_df,
    estimation_design = colnames(x$X1), time_step_hours = x$step,
    covariance = "sigma_g^2 D_g [(1-nugget) AR1(rho) + nugget I] D_g; independent cultures",
    precision = "D_g diagonal sqrt(1/wRNA + 1/wRibo); separate assay categorical-time voom",
    qualification = paste("F exact only for fixed correct precision/correlation and Gaussian errors.",
      "Estimated weights/rho/nugget uncertainty ignored; count log-bias and additive-AR covariance misspecification remain.",
      "Assay residual cross-covariance omitted from the working precision; gene residual scale is fitted."),
    moderate = moderate))
  ans
}

# Called after public input validation. An optimization warning does not cause
# silent replacement of finite p values; its convergence flag remains attached.
.cte_fit_gls <- function(data,
    method = c("weighted_ar", "unweighted_ar", "weighted_independent", "weighted_ar_ordinary"),
    period = 24, trend = TRUE) {
  method <- match.arg(method)
  args <- c(data[c("rna", "ribo", "meta", "off_rna", "off_ribo")],
            list(period = period, trend = trend, block_rep = TRUE))
  if (method == "unweighted_ar") {
    x <- do.call(.cte_ar_prepare, args)
    x$precision <- matrix(1, nrow(x$y), ncol(x$y))
    x$informative <- rowSums(data$rna) > 0 & rowSums(data$ribo) > 0
    x$precision_design <- character(0)
    x$precision_residual_df <- NA_integer_
  } else {
    x <- do.call(.cte_wg_prepare, args)
  }
  independent <- method == "weighted_independent"
  est <- if (independent) list(rho = 0, nugget = 1, converged = TRUE,
    genes = ncol(x$y), estimated = FALSE) else .cte_wg_estimate(x)
  moderate <- method != "weighted_ar_ordinary"
  ans <- .cte_wg_annotate(.cte_wg_test(x, est$rho, est$nugget, moderate),
                         x, est, method, moderate)
  if (method == "unweighted_ar")
    attr(ans, "diagnostics")$precision <-
      "All observation precisions fixed to 1; matched H1 covariance estimation"
  ans$fit_failed <- !is.finite(ans$p)
  ans$p[ans$fit_failed] <- 1
  attr(ans, "diagnostics")$calibrated <- FALSE
  attr(ans, "diagnostics")$engine <- "paired_logTE_GLS"
  ans
}

