#' @title Important principal components
#' @description Return number of important principal components
#'
#' @param pca principal component object.
#'
#' @return PCA object with only components that explain more variance than
#'   1 / number of variables.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
impPCs <- function(pca) {
  # number of important PCs = as many as account for expected variance
  imp.gt.exp <- pca$importance["Proportion of Variance", ] >= (1 / ncol(pca$importance))
  n.imp <- max(1, which.min(imp.gt.exp) - 1)
  
  # extract important PCs from object
  pca$rotation <- pca$rotation[, 1:n.imp]
  pca$x<- pca$x[1:n.imp, 1:n.imp]
  pca$importance <- pca$importance[, 1:n.imp]
  
  # remove unnecessary elements
  pca$sdev <- pca$center <- pca$scale <- NULL
  
  pca
}
