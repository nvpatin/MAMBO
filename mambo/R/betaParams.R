#' @title Fit beta distribution parameters
#' @description Fits shape paraemters of beta distribution to ASV occurrence data.
#'
#' @param x either a filename of a tab-separated file containing ASV occurrence
#'   data or a data frame of occurrence data. Occurrence data should have ASVs as 
#'   rows and samples as columns. 
#'
#' @return a three dimensional array where the first dimension are the samples,
#'   the second dimension are the ASVs, and the third dimension are the two beta 
#'   shape parameters.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
betaParams <- function(x) {
  # read data if filename is given
  if(is.character(x)) x <- utils::read.delim(x, row.names = 1)
  x <- as.matrix(x)

  if(length(rownames(x)) == 0) stop('count tables need ASV names for the rows')
  if(length(colnames(x)) == 0) stop('count tables need sample names for the columns')
  
  # by-sample (columns) coverage
  coverage <- colSums(x)

  # fit beta shape parameters (transpose matrix for recycling of coverage vector)
  tx <- t(x)
  result <- array(c(tx + 1, coverage - tx + 1), dim = c(dim(tx), 2))
  dimnames(result) <- list(colnames(x), rownames(x), c('shape1', 'shape2'))
  
  result
}
