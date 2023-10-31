#' @title Fit beta distribution parameters
#' @description Fits shape paraemters of beta distribution to ASV occurrence data.
#'
#' @param x either a filename of a tab-separated file containing ASV occurrence
#'   data or a matrix of occurrence data.
#'
#' @return an array with parameters of beta distribution for the occurrence of
#'   each ASV for each sample.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' # Random occurrence matrix
#'
#'
#' @export
#'
betaParams <- function(x) {
  # read data if filename is given
  if(is.character(x)) {
    x <- utils::read.delim(x, row.names = 1)
  }

  # by-sample (columns) coverage
  coverage <- colSums(x)

  # fit beta shape parameters (transpose matrix for recycling of coverage vector)
  tx <- t(x)
  result <- array(c(tx + 1, coverage - tx + 1), dim = c(dim(tx), 2))
  dimnames(result) <- list(
    colnames(x),
    rownames(x),
    c('shape1', 'shape2')
  )
  result
}
