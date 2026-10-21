#' Small synthetic paired example
#'
#' Generates 120 genes at ten four-hourly times in two cultures, with a known
#' illustrative TE effect in 20 genes. This is an API demonstration, not the
#' scientific benchmark simulator and not a power or FDR validation dataset.
#' The caller's random-number state is preserved.
#' @param seed Integer random seed.
#' @return A te_data object with simulation origin.
#' @export
example_te_data <- function(seed = 1L) {
  if (!is.numeric(seed) || length(seed)!=1L || !is.finite(seed) || seed < 0 || seed > .Machine$integer.max || seed != round(seed))
    stop("seed must be one nonnegative integer")
  had <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  if (had) old <- get(".Random.seed", envir=.GlobalEnv)
  on.exit(if(had) assign(".Random.seed",old,envir=.GlobalEnv) else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv))
  set.seed(seed)
  time <- rep(seq(12,48,4),2); batch <- rep(c("culture1","culture2"),each=10)
  mu <- outer(exp(stats::runif(120,log(30),log(500))), exp(.3*cos(2*pi*time/24)))
  effect <- matrix(0,120,20); effect[1:20,] <- rep(cos(2*pi*time/24),each=20)
  rna <- matrix(stats::rnbinom(length(mu),mu=as.vector(mu),size=20),120,20)
  ribo <- matrix(stats::rnbinom(length(mu),mu=as.vector(mu*2^effect),size=20),120,20)
  ids <- paste0("sample",seq_along(time)); gn <- sprintf("gene%03d",1:120)
  dimnames(rna) <- dimnames(ribo) <- list(gn,ids)
  te_data(rna,ribo,data.frame(sample_id=ids,time_hours=time,batch=batch),
    origin="simulation",dataset_id=paste0("example-",seed),calibration_context="illustrative_example_v1")
}
