.cte_test_gls_data <- function() {
  set.seed(73019)
  meta <- data.frame(time_hours = rep(seq(12, 48, 4), each = 2),
                     replicate = rep(1:2, 10), bio_sample_id = paste0("s", 1:20))
  base <- rep(c(20, 80, 320, 1280), each = 30)
  rnamp <- rep(c(0, 1, 2), 40)
  phase <- rep(c(0, 6, 12, 18), 30)
  rnwave <- cos(outer(phase, meta$time_hours, function(a,b) 2*pi*(b-a)/24))
  mu <- base * 2^(rnamp * rnwave)
  te <- rep(c(0, .8), each = 60) * sin(2*pi*meta$time_hours/24)[col(mu)]
  rna <- matrix(stats::rnbinom(120*20, mu = mu, size = 10), 120, 20)
  ribo <- matrix(stats::rnbinom(120*20, mu = mu*.6*2^te, size = 8), 120, 20)
  ribo[1,] <- 0; rna[2,] <- 0; rna[3,] <- ribo[3,] <- 0
  rownames(rna) <- rownames(ribo) <- paste0("g", 1:120)
  colnames(rna) <- colnames(ribo) <- meta$bio_sample_id
  list(rna = rna, ribo = ribo, meta = meta,
       off_rna = log(colSums(rna)), off_ribo = log(colSums(ribo)))
}

test_that("GLS returns stable probabilities and conservative whole-zero rows", {
  data <- .cte_test_gls_data()
  for (method in c("weighted_ar", "unweighted_ar", "weighted_independent", "weighted_ar_ordinary")) {
    fit <- .cte_fit_gls(data, method)
    expect_identical(fit$gene_id, rownames(data$rna))
    expect_true(all(is.finite(fit$p) & fit$p >= 0 & fit$p <= 1))
    expect_false(any(fit$fit_failed))
    expect_identical(fit$informative[1:3], rep(FALSE, 3))
    expect_identical(fit$p[1:3], rep(1, 3))
    expect_false(attr(fit, "diagnostics")$calibrated)
    expect_match(attr(fit, "diagnostics")$qualification, "uncertainty ignored")
  }
})

test_that("GLS retains the frozen small-simulation numerical regression", {
  # Captured from the frozen original, independently checked against the port.
  # A dependency algorithm change requires reviewing, not silently updating,
  # these expectations. This is an implementation test, not FDR validation.
  expected <- list(
    weighted_ar = c(.0330419684770636, .00330009522689469, .00397071750859112,
      .00712245019430158, .163022332682435, .205157432150656, .504727206757273,
      -.217842893008326, -.921648793221688, -1.14940061893488, -1.17901309138594),
    unweighted_ar = c(.072202255116438, .0126111748976932, .0233740404010072,
      .00196375301908736, .121421630478203, .202999163548848, .510476824608371,
      -.185573050982234, -.881221579030646, -1.16358256781215, -1.16959793980153),
    weighted_independent = c(.0253610508150011, .00264086599754833, .00269785143765074,
      .00472906258964859, .151768277565408, .205375941498253, .503480621779387,
      -.216215653175861, -.917663221347197, -1.12283456176486, -1.18289927014836),
    weighted_ar_ordinary = c(.125182181249401, .0461397613656976, .0495446435784123,
      .0135175231943783, .140703421120514, .205157432150656, .504727206757273,
      -.217842893008326, -.921648793221688, -1.14940061893488, -1.17901309138594))
  data <- .cte_test_gls_data()
  for (method in names(expected)) {
    a <- .cte_fit_gls(data, method)
    expect_equal(c(a$p[4:8], a$cos[4:6], a$sin[4:6]), expected[[method]], tolerance = 1e-5)
  }
})

test_that("GLS retains irregular-time and rank-deficiency checks", {
  data <- .cte_test_gls_data()
  bad <- data; bad$meta$time_hours[] <- 12
  expect_error(.cte_fit_gls(bad), "unique time points")
  bad <- data; bad$meta$time_hours <- rep(seq(0, 216, 24), each = 2)
  expect_error(.cte_fit_gls(bad), "full rank")
  bad <- data; bad$meta$replicate[] <- 1
  expect_error(.cte_fit_gls(bad), "unique time points")
  single <- data; keep <- seq(1, 19, 2)
  single$rna <- data$rna[,keep]; single$ribo <- data$ribo[,keep]
  single$meta <- data$meta[keep,]; single$off_rna <- data$off_rna[keep]
  single$off_ribo <- data$off_ribo[keep]
  expect_error(.cte_fit_gls(single, "weighted_independent"), "Precision design needs residual df")
  irregular <- data; irregular$meta$time_hours[3:4] <- 16.5
  fit <- .cte_fit_gls(irregular)
  expect_true(all(is.finite(fit$p)))
  expect_gte(attr(fit, "diagnostics")$rho, 0)
})

test_that("observation order does not change GLS estimates", {
  data <- .cte_test_gls_data(); perm <- c(seq(2,20,2), seq(1,19,2))
  other <- data
  other$rna <- data$rna[,perm]; other$ribo <- data$ribo[,perm]
  other$meta <- data$meta[perm,]; other$off_rna <- data$off_rna[perm]
  other$off_ribo <- data$off_ribo[perm]
  a <- .cte_fit_gls(data); b <- .cte_fit_gls(other)
  expect_equal(a$p, b$p, tolerance = 1e-5)
  expect_equal(a$cos, b$cos, tolerance = 1e-5)
  expect_equal(a$sin, b$sin, tolerance = 1e-5)
})

test_that("identifiable six-hour H1 schedules do not require H2 identification", {
  original <- .cte_test_gls_data()
  for (n_times in c(10, 4)) {
    keep <- seq_len(2*n_times)
    data <- original
    data$rna <- original$rna[,keep]; data$ribo <- original$ribo[,keep]
    data$meta <- original$meta[keep,]
    data$meta$time_hours <- rep(seq(0, by = 6, length.out = n_times), each = 2)
    data$off_rna <- original$off_rna[keep]; data$off_ribo <- original$off_ribo[keep]
    for (method in c("weighted_ar", "unweighted_ar")) {
      fit <- .cte_fit_gls(data, method)
      expect_equal(nrow(fit), nrow(data$rna))
      expect_true(all(is.finite(fit$p) & fit$p >= 0 & fit$p <= 1))
      expect_true(all(is.finite(fit$cos) & is.finite(fit$sin)))
    }
  }
})

test_that("ported GLS equals the frozen historical implementation when available", {
  # This additional development check skips in a distributed package, where
  # historical project sources are intentionally not installed dependencies.
  roots <- c(".", "..", "../..", "../../..")
  hit <- roots[file.exists(file.path(roots,
    "circadian_te_lowcount_v3_20261021/methods/weighted_gls.R"))]
  skip_if(length(hit) == 0, "Historical source is outside distributed package")
  ref <- new.env(parent = baseenv())
  # Package imports used unqualified by the frozen source.
  for (nm in c("model.matrix", "lm.fit", "optim", "pf")) ref[[nm]] <- get(nm, asNamespace("stats"))
  sys.source(file.path(hit[1], "circadian_te_optimization_v2_20261020/methods/ar_models.R"), envir = ref)
  sys.source(file.path(hit[1], "circadian_te_lowcount_v3_20261021/methods/weighted_gls.R"), envir = ref)
  data <- .cte_test_gls_data()
  original <- do.call(ref$fit_weighted_family, data)
  map <- c(weighted_ar = "weighted_ar_moderated", unweighted_ar = "unweighted_ar_h1_moderated",
           weighted_independent = "weighted_independent", weighted_ar_ordinary = "weighted_ar")
  for (method in names(map)) {
    fit <- .cte_fit_gls(data, method); old <- original[[map[[method]]]]
    for (col in c("p", "cos", "sin", "amplitude_log2", "phase_hours"))
      expect_equal(fit[[col]], old[[col]], tolerance = 1e-10, info = paste(method, col))
  }
})
