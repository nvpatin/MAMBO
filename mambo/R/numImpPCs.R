#' @title Number of important principal components
#' @description Compute number of important principal components
#'
#' @param pca principal component object.
#'
#' @return number of components that explain more variance than
#'   1 / number of variables.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' # Random occurrence matrix
#'
#' @export
#'
numImpPCs <- function(pca) {
  # number of important PCs = as many as account for expected variance
  imp.gt.exp <- pca$importance["Proportion of Variance", ] >= (1 / ncol(pca$importance))
  max(1, which.min(imp.gt.exp) - 1)
}
