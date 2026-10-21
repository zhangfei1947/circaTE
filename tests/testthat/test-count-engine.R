.count_fixture <- function(ng = 160L) {
  set.seed(846101)
  meta <- data.frame(time_hours = rep(seq(12, 48, 4), each = 2),
    replicate = rep(1:2, 10), bio_sample_id = paste0("s", seq_len(20)))
  mu <- exp(stats::runif(ng, log(25), log(800))) *
    2^cos(outer(stats::runif(ng, 0, 2 * pi), 2 * pi * meta$time_hours / 24, "-"))
  r <- matrix(stats::rnbinom(ng * 20, mu = mu, size = 15), ng, 20)
  b <- matrix(stats::rnbinom(ng * 20, mu = .7 * mu, size = 10), ng, 20)
  rownames(r) <- rownames(b) <- paste0("g", seq_len(ng))
  colnames(r) <- colnames(b) <- meta$bio_sample_id
  list(rna = r, ribo = b, meta = meta,
       off_rna = rep(log(1e6), 20), off_ribo = rep(log(1e6), 20))
}

test_that("count engines distinguish uninformative genes from fitting failures", {
  skip_if_not_installed("edgeR", minimum_version = "4.4.0")
  d <- .count_fixture()
  d$rna[1, ] <- 0; d$ribo[1, ] <- 0
  d$rna[2, ] <- 0
  d$ribo[3, ] <- 0
  methods <- "edger"
  if (requireNamespace("DESeq2", quietly = TRUE)) methods <- c(methods, "deseq2")
  for (m in methods) {
    z <- .cte_fit_count(d, method = m)
    expect_identical(z$gene_id, rownames(d$rna))
    expect_equal(z$p[1:3], rep(1, 3))
    expect_false(any(z$informative[1:3]))
    expect_false(any(z$fit_failed[1:3]))
    expect_true(all(is.na(z$amplitude_log2[1:3])))
    expect_identical(z$status[1:3], rep("uninformative", 3))
    expect_true(all(is.finite(z$p) & z$p >= 0 & z$p <= 1))
    expect_equal(attr(z, "diagnostics")$genes_uninformative, 3)
    # Removing uninformed genes before fitting should not change the remaining fits.
    dd <- d; dd$rna <- d$rna[-(1:3), ]; dd$ribo <- d$ribo[-(1:3), ]
    zz <- .cte_fit_count(dd, method = m)
    expect_equal(z$p[-(1:3)], zz$p, tolerance = 1e-12)
  }
})

test_that("all-zero input is retained without pretending optimizer failures", {
  skip_if_not_installed("edgeR", minimum_version = "4.4.0")
  d <- .count_fixture(4L); d$rna[,] <- 0; d$ribo[,] <- 0
  z <- .cte_fit_count(d, "edger")
  expect_equal(z$p, rep(1, 4))
  expect_false(any(z$informative | z$fit_failed))
  expect_true(attr(z, "diagnostics")$no_informative_genes)
})

test_that("internal count errors are explicit", {
  skip_if_not_installed("edgeR", minimum_version = "4.4.0")
  d <- .count_fixture(10L)
  expect_error(.cte_fit_count(d, "unknown"), "arg")
  expect_error(.cte_fit_count(d, "edger", period = 0), "period")
  bad <- d; bad$rna[1, 1] <- -.5
  expect_error(.cte_fit_count(bad, "edger"), "integer")
  bad <- d; bad$off_rna[1] <- Inf
  expect_error(.cte_fit_count(bad, "edger"), "offset")
  bad <- d; bad$meta$time_hours[] <- 12
  expect_error(.cte_fit_count(bad, "edger"), "rank deficient")
  bad <- d; bad$rna[-1, ] <- 0
  expect_error(.cte_fit_count(bad, "edger"), "two informative")
})

test_that("packaged count engines numerically reproduce frozen 720-gene engines", {
  skip_if_not_installed("edgeR", minimum_version = "4.4.0")
  # Optional checkout regression: package consumers do not need historical data.
  ancestors <- vapply(0:4, function(k) paste(c(".", rep("..", k)), collapse = "/"), character(1))
  root <- ancestors[file.exists(file.path(ancestors,
    "circadian_te_modeling_20261020/methods/count_models.R"))]
  skip_if(length(root) == 0L, "Frozen development checkout not present")
  root <- root[1]
  input <- file.path(root, "circadian_te_lowcount_v3_20261021/inputs/simulations",
                     "development__01__rho0_notch4.rds")
  skip_if_not(file.exists(input), "Frozen 720-gene development fixture not present")
  g <- readRDS(input)
  expect_equal(nrow(g$rna), 720)
  old <- new.env(parent = globalenv())
  sys.source(file.path(root, "circadian_te_modeling_20261020/methods/count_models.R"), old)
  sys.source(file.path(root, "circadian_te_modeling_20261020/methods/deseq_model.R"), old)
  sys.source(file.path(root, "circadian_te_optimization_v2_20261020/methods/deseq_calibration.R"), old)
  args <- g[c("rna", "ribo", "meta", "off_rna", "off_ribo")]
  cols <- c("p", "cos", "sin", "amplitude_log2", "phase_hours")
  reference <- do.call(old$fit_edger_paired, args)
  now <- .cte_fit_count(g, "edger")
  expect_equal(now[, cols], reference[, cols], tolerance = 1e-10)
  expect_false(any(now$fit_failed | !now$informative))
  expect_identical(attr(now, "diagnostics")$legacy, FALSE)
  if (requireNamespace("DESeq2", quietly = TRUE)) {
    reference <- do.call(old$fit_deseq_paired, args)
    now <- .cte_fit_count(g, "deseq2")
    expect_equal(now[, cols], reference[, cols], tolerance = 1e-10)
    finite <- .cte_fit_count(g, "deseq2_finite")
    old_finite <- old$deseq_finite_table(reference)
    expect_equal(finite[, cols], old_finite[, cols], tolerance = 1e-10)
    expect_true(all(finite$p >= now$p - 1e-14))
    expect_false(attr(finite, "diagnostics")$calibration_exact)
  }
})
