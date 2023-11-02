rm(list = ls())
library(mambo)

results <- readRDS('mambo.20231031_112830.rds')

pca.l <- extractPCA(results)$loadings


apply(pca$rotation[, 1:pca$num.pcs], 2, function(x) {
  outliers <- x[isOutlier(x)]
  list(
    pos = sort(outliers[outliers > 0], decreasing = TRUE),
    neg = sort(outliers[outliers < 0]))
})