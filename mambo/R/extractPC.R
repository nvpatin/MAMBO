#' @title Extract principal component loadings and scores
#' @description Separate and sort positive and negative components
#'
#' @param pca.list list of PCA results.
#' @param pc principal component to summarize (1:n).
#' @param locus name or number of locus to extract.
#'
#' @return a list with component loadings (\code{$loadings}) and scores
#'   (\code{$scores}).
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
extractPC <- function(pca.list, pc, locus) {
  loadings <- sapply(pca.list, function(x) x$pca[[locus]]$rotation[, pc])
  for(i in 1:nrow(loadings)) {
    switch.sign <- sign(loadings[i, ]) != sign(loadings[i, 1])
    loadings[i, switch.sign] <- loadings[i, switch.sign] * -1
  }

  scores <- sapply(pca.list, function(x) x$pca[[locus]]$x[, pc])
  for(i in 1:nrow(scores)) {
    switch.sign <- sign(scores[i, ]) != sign(scores[i, 1])
    scores[i, switch.sign] <- scores[i, switch.sign] * -1
  }

  list(loadings = loadings, scores = scores)
}
