#' @title Fit beta distribution parameters
#' @description Fits shape paraemters of beta distribution to ASV occurrence data.
#'
#' @param x either a filename of a tab-separated file containing ASV occurrence
#'   data or a matrix of occurrence data.
#'
#' @return a three dimensional array where the first dimension are the samples,
#'   the second dimension are the ASVs, and the third dimension are the two beta 
#'   shape parameters
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' # use built-in 16S data
#' beta.16s <- betaParams(fl16s)
#' str(beta.16s)
#' 
#' # both shape parameters for first 5 samples and first 3 ASVs
#' beta.16s[1:5, 1:3, ]
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
