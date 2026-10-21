# Independent public-API regression checks from the count-engine review.
.review_count_fit <- local({
  cache <- new.env(parent = emptyenv())
  function(seed = 9101L) {
    key <- as.character(seed)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      d <- example_te_data(seed)
      # These fits exercise validation only; no statistical all-null claim is made.
      assign(key, fit_te(d, "edger"), envir = cache)
    }
    get(key, envir = cache, inherits = FALSE)
  }
})

test_that("public fitting revalidates mutable input and keeps single-assay zeros", {
  d <- example_te_data(9103L)
  bad <- d; bad$metadata$sample_id[1] <- "wrong-id"
  expect_error(fit_te(bad, "edger"), "sample_id")
  bad <- d; bad$rna[1, 1] <- -1
  expect_error(fit_te(bad, "edger"), "nonnegative")
  d$rna[1, ] <- 0; d$ribo[2, ] <- 0
  fit <- fit_te(d, "edger")
  expect_equal(fit$results$p[1:2], c(1, 1))
  expect_equal(fit$results$q[1:2], c(1, 1))
  expect_false(any(fit$results$informative[1:2]))
  expect_false(any(fit$results$fit_failed[1:2]))
})

test_that("calibration rejects edited result alignment and edited designs", {
  a <- .review_count_fit(9101L); b <- .review_count_fit(9102L)
  bad <- b; bad$results <- bad$results[1:2, ]
  expect_error(build_te_calibrator(list(a, bad), c("x", "x"), all_null = TRUE))
  bad <- b; bad$results$gene_id[1] <- "not-the-input-gene"
  expect_error(build_te_calibrator(list(a, bad), c("x", "x"), all_null = TRUE))
  bad <- b; bad$data$metadata$time_hours[1] <- bad$data$metadata$time_hours[1] + 1
  expect_error(build_te_calibrator(list(a, bad), c("x", "x"), all_null = TRUE))
})

test_that("renaming a dataset does not make reused calibration observations independent", {
  a <- .review_count_fit(9101L)
  duplicate <- a; duplicate$data$dataset_id <- "same-counts-renamed"
  expect_error(build_te_calibrator(list(a, duplicate), c("x", "x"), all_null = TRUE))
  b <- .review_count_fit(9102L)
  ca <- build_te_calibrator(list(a, b), c("x", "x"), all_null = TRUE)
  expect_error(calibrate_te(duplicate, ca))
})

test_that("QC run lengths use within-culture sampling intervals", {
  d <- example_te_data(9104L)
  md <- d$metadata
  md$time_hours[md$batch == "culture2"] <- md$time_hours[md$batch == "culture2"] + 2
  r <- d$rna; b <- d$ribo; r[1, ] <- 1; b[1, ] <- 1
  z <- te_data(r, b, md, offsets = list(rna = d$off_rna, ribo = d$off_ribo))
  qc <- te_qc(z)$genes
  # Both cultures have10consecutive4h observations despite an offset between them.
  expect_equal(qc$longest_low_rna_run[1], 10L)
  expect_equal(qc$longest_low_ribo_run[1], 10L)
})
