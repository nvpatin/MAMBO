#' @title Random relative percent occurrence
#' @description Random sample of percent occurrence from beta parameters.
#'
#' @param beta.params matrix of beta parameters from \link{betaParams}.
#' @param num.cores number of cores to use for sampling relative percent 
#'    occurrence. Defaults to value reported by \link[parallel]{detectCores} 
#'    minus 1.
#'
#' @return matrix of random percent occurrences. Rows are ASVs and columns 
#'   are samples.
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#' 
ranRelPct <- function(beta.params, num.cores = parallel::detectCores(logical = FALSE) - 1) {
  if(is.null(beta.params)) return(NULL)
  param.dim <- dim(beta.params)[1:2]
  # draw random sample from beta distribution with shape parameters 'p'
  pct <- parallel::mclapply(1:prod(param.dim), function(i) {
    inds <- arrayInd(i, param.dim)
    p <- beta.params[inds[1], inds[2], ]
    stats::rbeta(1, p[1], p[2])
  }, mc.cores = num.cores) |>
    unlist() |>
    array(dim = param.dim)
  
  # normalize random percents to unity and return matrix
  # (transposed back to original dimensions)
  pct <- t(pct / rowSums(pct))
  rownames(pct) <- dimnames(beta.params)[[2]]
  colnames(pct) <- dimnames(beta.params)[[1]]
  pct
}
