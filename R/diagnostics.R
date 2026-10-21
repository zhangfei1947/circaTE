#' Quantify low-count observations and contiguous troughs
#' @param data A te_data or te_fit object.
#' @param threshold Counts below this nonnegative threshold are considered low.
#' @return A list with gene and sample QC data frames. Runs are measured within
#'   each culture in time order; gaps larger than the smallest within-culture time step
#'   break runs, so missing observations are not counted as measured troughs.
#' @export
te_qc <- function(data, threshold = 10) {
  if (inherits(data, "te_fit")) data <- data$data
  if (!inherits(data, "te_data")) stop("data must be te_data or te_fit")
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) || threshold < 0)
    stop("threshold must be one nonnegative number")
  longest <- function(x) {
    out <- integer(nrow(x))
    for (b in unique(data$meta$replicate)) {
      j <- which(data$meta$replicate == b); j <- j[order(data$meta$time_hours[j])]
      step <- if (length(j) > 1L) min(diff(data$meta$time_hours[j])) else Inf
      run <- integer(nrow(x)); mx <- run
      for (k in seq_along(j)) {
        if (k > 1L && data$meta$time_hours[j[k]] - data$meta$time_hours[j[k-1L]] > step + 1e-7)
          run[] <- 0L
        run <- ifelse(x[, j[k]] < threshold, run + 1L, 0L)
        mx <- pmax(mx, run)
      }
      out <- pmax(out, mx)
    }
    out
  }
  g <- data.frame(gene_id = rownames(data$rna),
    mean_rna_count = rowMeans(data$rna), mean_ribo_count = rowMeans(data$ribo),
    min_rna_count = apply(data$rna, 1, min), min_ribo_count = apply(data$ribo, 1, min),
    n_low_rna = rowSums(data$rna < threshold), n_low_ribo = rowSums(data$ribo < threshold),
    longest_low_rna_run = longest(data$rna), longest_low_ribo_run = longest(data$ribo),
    informative = rowSums(data$rna) > 0 & rowSums(data$ribo) > 0)
  s <- data.frame(sample_id = data$metadata$sample_id, time_hours = data$meta$time_hours,
    batch = data$meta$replicate, rna_library_size = colSums(data$rna), ribo_library_size = colSums(data$ribo),
    rna_effective_library_size = exp(data$off_rna), ribo_effective_library_size = exp(data$off_ribo),
    n_low_rna_genes = colSums(data$rna < threshold), n_low_ribo_genes = colSums(data$ribo < threshold))
  list(genes = g, samples = s, threshold = threshold,
       definition = "Consecutive observed low counts within culture; gaps larger than the within-culture minimum step break runs (operational count, not duration).")
}

#' Plot paired gene trajectories
#' @param x A te_data or te_fit object.
#' @param gene_id One gene identifier.
#' @param pseudocount Positive pseudocount for the descriptive log ratio only.
#' @param ... Reserved for future use.
#' @return Invisibly a sample-level data frame. Draws normalized RNA, Ribo and
#'   log2 TE trajectories without claiming a significant fitted rhythm.
#' @export
plot_te_gene <- function(x, gene_id, pseudocount = 0.5, ...) {
  data <- if (inherits(x, "te_fit")) x$data else x
  if (!inherits(data, "te_data")) stop("x must be te_data or te_fit")
  if (length(gene_id) != 1L || is.na(gene_id) || !gene_id %in% rownames(data$rna)) stop("Unknown gene_id")
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L || !is.finite(pseudocount) || pseudocount <= 0)
    stop("pseudocount must be one positive number")
  i <- match(gene_id, rownames(data$rna))
  r <- data$rna[i, ]; b <- data$ribo[i, ]; lr <- exp(data$off_rna); lb <- exp(data$off_ribo)
  out <- data.frame(sample_id = data$metadata$sample_id, time_hours = data$meta$time_hours,
    batch = data$meta$replicate, rna_count = r, ribo_count = b,
    rna_cpm = 1e6*r/lr, ribo_cpm = 1e6*b/lb,
    log2_te = log2((b+pseudocount)/(lb+2*pseudocount)) - log2((r+pseudocount)/(lr+2*pseudocount)))
  old <- graphics::par(mfrow = c(1,3), mar=c(4,4,3,1)); on.exit(graphics::par(old))
  batches <- unique(out$batch); cols <- grDevices::hcl.colors(length(batches), "Dark 3")
  for (nm in c("rna_cpm", "ribo_cpm", "log2_te")) {
    graphics::plot(out$time_hours, out[[nm]], type="n", xlab="Time (hours)",
      ylab=nm, main=paste(gene_id, nm))
    for (k in seq_along(batches)) {
      z <- out[out$batch == batches[k], ]; z <- z[order(z$time_hours), ]
      graphics::lines(z$time_hours, z[[nm]], type="b", col=cols[k], pch=15+k)
    }
    if (nm == "rna_cpm") graphics::legend("topright", legend=batches, col=cols, pch=16+seq_along(batches)-1L, bty="n", cex=.75)
  }
  invisible(out)
}

#' Export a reproducible fit bundle
#' @param fit A te_fit object.
#' @param path A new output directory; existing nonempty directories are rejected.
#' @return Invisibly the normalized output path. Writes results/QC TSV, complete
#'   fit RDS, manifest JSON and session information. Files are not published.
#' @export
write_te_results <- function(fit, path) {
  if (!inherits(fit, "te_fit")) stop("fit must be te_fit")
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) stop("path must be one directory")
  if (file.exists(path) && (!dir.exists(path) || length(list.files(path, all.files=TRUE, no..=TRUE))))
    stop("Output path must be absent or an empty directory")
  if (!dir.exists(path) && !dir.create(path, recursive=TRUE)) stop("Could not create output directory")
  qc <- te_qc(fit)
  for (nm in c("results", "gene_qc", "sample_qc")) {
    tab <- switch(nm, results=fit$results, gene_qc=qc$genes, sample_qc=qc$samples)
    utils::write.table(tab, file.path(path,paste0(nm,".tsv")), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  }
  saveRDS(fit, file.path(path,"fit.rds"))
  jsonlite::write_json(list(spec=fit$spec, versions=fit$versions,
    inferential_status=fit$inferential_status, notes=fit$notes,
    diagnostics=.cte_plain_json(fit$diagnostics), calibration=fit$calibration,
    origin=fit$data$origin, dataset_id=fit$data$dataset_id,
    written_utc=format(Sys.time(),tz="UTC",usetz=TRUE)),
    file.path(path,"manifest.json"), pretty=TRUE, auto_unbox=TRUE, null="null", na="null")
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(path,"sessionInfo.txt"))
  writeLines(c("circaTE development analysis", "", fit$notes,
    "", "q values use all attempted genes in this single fit; q<0.05 is nominal.",
    "Analyses across groups, methods or periods do not share automatic multiple-testing correction.",
    "Source counts are retained in fit.rds. Low observations were not removed."), file.path(path,"README.txt"))
  invisible(normalizePath(path))
}

.cte_plain_json <- function(x) {
  if (is.list(x)) return(lapply(x, .cte_plain_json))
  if (is.atomic(x) && is.object(x)) return(unclass(x))
  x
}
