#!/usr/bin/env Rscript
# The installed package is required; no source-tree or working-directory assumptions.
args <- commandArgs(trailingOnly=TRUE)
help <- paste(
  "circaTE 0.1.0",
  "Usage: Rscript circaTE.R --rna RNA.tsv --ribo Ribo.tsv --metadata samples.tsv --out results",
  "Optional: --method edger|weighted_ar|unweighted_ar|weighted_independent|weighted_ar_ordinary|deseq2|deseq2_finite",
  "          --period 24 --normalization TMM|library_size --no-trend",
  "Counts: gene_id first column, biological sample IDs as remaining column names.",
  "Metadata: sample_id, time_hours, batch (one group per run).",
  "Output directory must be new or empty. All reported inference is exploratory.",sep="\n")
if (!length(args) || "--help" %in% args) {cat(help,"\n");quit(status=0)}
main <- function(args) {
  if (!requireNamespace("circaTE",quietly=TRUE)) stop("Install circaTE first")
  opt <- list(method="edger",period="24",normalization="TMM",trend=TRUE)
  keys <- c("rna","ribo","metadata","out","method","period","normalization")
  i <- 1L
  while (i<=length(args)) {
    a <- args[i]
    if (a=="--no-trend") {opt$trend <- FALSE;i <- i+1L;next}
    k <- sub("^--","",a)
    if (!startsWith(a,"--") || !k %in% keys || i==length(args) || startsWith(args[i+1L],"--"))
      stop("Unknown or incomplete option: ",a)
    opt[[k]] <- args[i+1L];i <- i+2L
  }
  if (!all(c("rna","ribo","metadata","out") %in% names(opt))) stop("Required options missing; use --help")
  d <- circaTE::read_te_data(opt$rna,opt$ribo,opt$metadata,normalization=opt$normalization)
  f <- circaTE::fit_te(d,method=opt$method,period=suppressWarnings(as.numeric(opt$period)),trend=opt$trend)
  circaTE::write_te_results(f,opt$out)
  print(f);cat("Saved:",normalizePath(opt$out),"\n")
}
tryCatch(main(args),error=function(e) {message("circaTE error: ",conditionMessage(e));quit(status=1)})
