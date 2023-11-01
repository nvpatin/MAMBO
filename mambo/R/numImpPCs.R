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
#' # use built-in 16S data
#' beta.16s <- betaParams(fl16s)
#' pca.16s <- ranPCA(beta.16s)
#' 
#' i <- numImpPCs(pca.16s)
#' 
#' # Number of variables
#' ncol(pca.16s$importance)
#' 
#' # Expected proportion of variance
#' 1 / ncol(pca.16s$importance)
#' 
#' # Number of important components
#' i
#' 
#' pca.16s$importance[, 1:i]
#'
#' @export
#'
numImpPCs <- function(pca) {
  # number of important PCs = as many as account for expected variance
  imp.gt.exp <- pca$importance["Proportion of Variance", ] >= (1 / ncol(pca$importance))
  max(1, which.min(imp.gt.exp) - 1)
}
