#' @title Random relative percent occurrence
#' @description Random sample of percent occurrence from beta parameters.
#'
#' @param beta.params matrix of beta parameters from \link{betaParams}.
#'
#' @return matrix of random percent occurrences. Rows are ASVs and columns 
#'   are samples.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' # use built-in 16S data
#' beta.16s <- betaParams(fl16s)
#' str(beta.16s)
#' 
#' rel.pct <- ranRelPct(beta.16s)
#' rel.pct[1:5, 1:3]
#'
#' @export
#'
ranRelPct <- function(beta.params) {
  # draw random sample from beta distribution with shape parameters 'p'
  pct <- apply(beta.params, c(1, 2), function(p) stats::rbeta(1, p[1], p[2]))

  # normalize random percents to unity and return matrix
  # (transposed back to original dimensions)
  pct <- t(pct / rowSums(pct))
  rownames(pct) <- dimnames(beta.params)[[2]]
  colnames(pct) <- dimnames(beta.params)[[1]]
  pct
}
