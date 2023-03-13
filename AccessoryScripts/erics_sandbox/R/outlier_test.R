rm(list = ls())

data(iris)

iris.pc <- princomp(iris[, 1:4])

loadings <- iris.pc$loadings

# function to return vector of logicals identifying if each 
# predictor is an outlier
isOutlier <- function(x) {
  quarts <- quantile(x, c(0.25, 0.75))
  iqr <- diff(quarts)
  lower <- quarts[1] - 1.5 * iqr
  upper <- quarts[2] + 1.5 * iqr
  x <= lower | x >= upper
}

# results for first component
out <- isOutlier(loadings[, 1])
# names of outlier predictors for first component
rownames(loadings)[out]

# return list of outlier names for all components
outliers <- apply(loadings, 2, function(x) {
  names(x)[isOutlier(x)]
}, simplify = FALSE)
outliers
