#! /usr/local/bin/Rscript

#' Draws a random relative percent occurrence matrix based on beta
#'   distribution fit to binomial
ranRelPct <- function(x, num.cores = NULL) {
  # by-sample (columns) coverage
  coverage <- colSums(x)

  # fit beta shape parameters (transpose matrix for recycling of coverage vector)
  tx <- t(x)
  beta.params <- array(c(tx + 1, coverage - tx + 1), dim = c(dim(tx), 2))

  # draw one random sample from beta distribution with shape parameters 'p'
  sample.beta <- function(p) rbeta(1, p[1], p[2])

  # draw matrix of random percent occurrence
  # use parallel processing if available
  if(is.null(num.cores)) num.cores <- parallel::detectCores() - 1
  pct <- if(num.cores > 1) {
    # spread draws across cores
    cl <- if(.Platform$OS.type == "windows") {
      parallel::makePSOCKcluster(num.cores)
    } else {
      parallel::makeForkCluster(num.cores)
    }
    pct <- parallel::parApply(cl, beta.params, c(1, 2), sample.beta)
    parallel::stopCluster(cl)
    pct
  } else {
    # single core mode
    apply(beta.params, c(1, 2), sample.beta)
  }

  # normalize random percents to unity and return matrix
  # (transposed back to original dimensions)
  t(pct / rowSums(pct))
}
