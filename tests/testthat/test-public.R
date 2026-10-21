test_that("pair alignment and malformed counts fail closed", {
  d <- example_te_data()
  expect_error(te_data(d$rna,d$ribo[,20:1],d$metadata),"order")
  m <- d$metadata; m$sample_id <- rev(m$sample_id)
  expect_error(te_data(d$rna,d$ribo,m),"order")
  for(v in c(-1,NA,Inf,.25)) {
    r <- d$rna;r[1,1] <- v
    expect_error(te_data(r,d$ribo,d$metadata),"finite nonnegative integers")
  }
  r <- d$rna;r[,1] <- 0
  expect_error(te_data(r,d$ribo,d$metadata),"empty library")
  m <- d$metadata;m$group <- rep(c("A","B"),each=10)
  expect_error(te_data(d$rna,d$ribo,m),"one biological group")
  expect_error(te_data(d$rna,d$ribo,d$metadata,offsets=list(rna=rep(1e8,20),ribo=d$off_ribo)),"log effective")
})

test_that("normalization is explicit and QC runs break at gaps", {
  d <- example_te_data(); e <- te_data(d$rna,d$ribo,d$metadata,offsets=list(rna=d$off_rna,ribo=d$off_ribo))
  expect_equal(d$off_rna,e$off_rna)
  q <- te_qc(d); expect_equal(nrow(q$genes),120);expect_equal(nrow(q$samples),20)
  d$rna[1,] <- 20; d$rna[1,1:3] <- 1
  expect_equal(te_qc(d)$genes$longest_low_rna_run[1],3)
  j <- setdiff(1:20,2)
  e <- te_data(d$rna[,j],d$ribo[,j],d$metadata[j,],offsets=list(rna=d$off_rna[j],ribo=d$off_ribo[j]))
  expect_equal(te_qc(e)$genes$longest_low_rna_run[1],1)
})

test_that("public fit preserves family and rejects unsupported designs", {
  d <- example_te_data();d$rna[1,] <- 0;d$ribo[1,] <- 0
  f <- fit_te(d)
  expect_s3_class(f,"te_fit");expect_equal(nrow(f$results),120)
  expect_equal(f$results$p[1],1);expect_equal(f$results$status[1],"uninformative")
  expect_equal(f$results$q,p.adjust(f$results$p,"BH"))
  expect_error(fit_te(d,period=0),"positive")
  expect_error(fit_te(d,period=c(22,24)),"positive")
  expect_error(fit_te(d,trend=NA),"TRUE or FALSE")
  j <- 1:10;one <- te_data(d$rna[,j],d$ribo[,j],d$metadata[j,])
  expect_error(fit_te(one,method="weighted_ar"),"two independent")
  expect_error(fit_te(one,method="weighted_independent"),"residual|replicat|Precision")
  out <- tempfile();write_te_results(f,out)
  expect_true(file.exists(file.path(out,"manifest.json")))
  expect_equal(readRDS(file.path(out,"fit.rds"))$results,f$results)
  expect_error(write_te_results(f,out),"empty")
})

test_that("calibration rejects leakage and real data and respects inclusive ties", {
  d <- example_te_data();f <- fit_te(d)
  a <- fit_te(example_te_data(2)); b <- fit_te(example_te_data(3))
  a$data$dataset_id <- "null-a";b$data$dataset_id <- "null-b"
  # Unit fixture for mapping arithmetic, not a scientific null simulation.
  a$results$p <- rep(c(0.01,.1,.1,1),30);b$results$p <- rep(c(.02,.2,.2,1),30)
  expect_error(build_te_calibrator(list(a,b),c("A","B")),"attestation")
  cal <- build_te_calibrator(list(a,b),c("A","B"),all_null=TRUE)
  expect_error(calibrate_te(a,cal),"overlaps")
  real <- f;real$data$origin <- "observed"
  expect_error(calibrate_te(real,cal),"simulation")
  f$results$p <- rep(c(0,.1,1),40)
  z <- calibrate_te(f,cal,"max")
  expected <- vapply(f$results$p,function(p) max(vapply(cal$strata,function(x) (1+sum(x<=p))/(length(x)+1),numeric(1))),numeric(1))
  expect_equal(z$results$p,expected);expect_equal(z$results$q,p.adjust(expected,"BH"))
  expect_error(calibrate_te(z,cal),"Already")
  dup <- a; dup$data$dataset_id <- "renamed"
  expect_error(calibrate_te(dup,cal),"counts overlap")
  short <- a; short$results <- short$results[1:2,]
  expect_error(build_te_calibrator(list(short,b),c("A","B"),TRUE),"complete aligned")
  f$spec$period <- 22
  expect_error(calibrate_te(f,cal),"does not match")
})

test_that("example generator preserves caller RNG and TSV roundtrip aligns", {
  set.seed(42);a <- .Random.seed;d <- example_te_data();expect_identical(.Random.seed,a)
  tmp <- tempfile();dir.create(tmp)
  for(nm in c("rna","ribo")) write.table(data.frame(gene_id=rownames(d[[nm]]),d[[nm]],check.names=FALSE),file.path(tmp,paste0(nm,".tsv")),sep="\t",row.names=FALSE,quote=FALSE)
  write.table(d$metadata,file.path(tmp,"metadata.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  x <- read_te_data(file.path(tmp,"rna.tsv"),file.path(tmp,"ribo.tsv"),file.path(tmp,"metadata.tsv"))
  expect_equal(x$rna,d$rna);expect_equal(x$ribo,d$ribo)
})


test_that("floating-point sine remnants do not rescue an aliased H1 design", {
  d <- example_te_data()
  d$metadata$time_hours <- rep(seq(0,108,12),2)
  expect_error(fit_te(d,"edger",period=24),"aliased")
  expect_error(fit_te(d,"weighted_ar",period=24),"aliased")
})
