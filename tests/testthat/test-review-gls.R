test_that("public GLS dispatch preserves offsets and the internal model", {
  d <- example_te_data(11)
  for (method in c("weighted_ar", "unweighted_ar", "weighted_independent", "weighted_ar_ordinary")) {
    public <- fit_te(d, method)
    internal <- .cte_fit_gls(d, method)
    expect_equal(public$results$p, internal$p, tolerance = 1e-12)
    expect_equal(public$results$q, stats::p.adjust(internal$p, "BH"))
    expect_identical(public$data$normalization, d$normalization)
    expect_identical(public$data$off_rna, d$off_rna)
    expect_identical(public$data$off_ribo, d$off_ribo)
    expect_false(public$diagnostics$calibrated)
  }
})

test_that("public refit rebuilds internal metadata and rejects broken alignment", {
  d <- example_te_data(12)
  broken <- d; broken$metadata$sample_id <- rev(broken$metadata$sample_id)
  expect_error(fit_te(broken, "weighted_ar"), "sample_id order")
  stale <- d; stale$meta$time_hours[] <- 0
  expect_equal(fit_te(stale, "weighted_ar")$results$p,
               fit_te(d, "weighted_ar")$results$p, tolerance = 1e-12)
})

test_that("descriptive default logTE agrees with the GLS response", {
  d <- example_te_data(13)
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  out <- tryCatch(plot_te_gene(d, "gene001"), finally = grDevices::dev.off())
  x <- .cte_ar_prepare(d$rna,d$ribo,d$meta,d$off_rna,d$off_ribo,24,TRUE,TRUE)
  expect_equal(out$log2_te, unname(x$y[,1]), tolerance = 1e-12)
  unlink(file)
})
