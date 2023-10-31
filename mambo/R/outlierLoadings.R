#' @title Outlier principal component loadings
#' @description Return loadings of outlier ASVs
#'
#' @param pca principal component object.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
outlierLoadings <- function(pca) {
  apply(pca$rotation[, 1:pca$num.pcs], 2, function(x) {
    outliers <- x[isOutlier(x)]
    list(
      pos = sort(outliers[outliers > 0], decreasing = TRUE),
      neg = sort(outliers[outliers < 0]))
  })
}
