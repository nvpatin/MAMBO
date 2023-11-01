#' @title Extract principal component loadings and scores
#' @description Extract principal component loadings and scores.
#'
#' @param results output of a \code{mambo} run.
#'
#' @return a list with a data frame of component loadings (\code{$loadings}) and scores
#'   (\code{$scores}) for each locus.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
extractPCA <- function(results) {
  pca.list <- results$reps |> 
    lapply(function(x) x$pca) |> 
    purrr::list_transpose() |> 
    lapply(purrr::list_transpose)
  
  loadings <- lapply(pca.list, function(x) {
    res <- abind::abind(x$rotation, along = 3) |> 
      apply(2, function(x) {
        to.switch <- sign(x[1, ]) != sign(x[1, 1])
        x[, to.switch] <- x[, to.switch] * -1
        x
      }, simplify = FALSE) |> 
      abind::abind(along = 3) |> 
      aperm(c(1, 3, 2))
    dimnames(res)[[3]] <- 1:dim(res)[3]
    as.data.frame.table(res) |> 
      stats::setNames(c('sample', 'pc', 'rep', 'loading')) |> 
      dplyr::mutate(
        sample = as.character(sample),
        pc = as.numeric(gsub('PC', '', pc)),
        rep = as.numeric(rep)
      )
  })
  
  scores <- lapply(pca.list, function(x) {
    res <- abind::abind(x$x, along = 3) |> 
      apply(2, function(x) {
        to.switch <- sign(x[1, ]) != sign(x[1, 1])
        x[, to.switch] <- x[, to.switch] * -1
        x
      }, simplify = FALSE) |> 
      abind::abind(along = 3) |> 
      aperm(c(1, 3, 2))
    dimnames(res)[[3]] <- 1:dim(res)[3]
    as.data.frame.table(res) |> 
      stats::setNames(c('sample', 'pc', 'rep', 'score')) |> 
      dplyr::mutate(
        sample = as.character(sample),
        pc = as.numeric(gsub('PC', '', pc)),
        rep = as.numeric(rep)
      )
  })
  
  list(loadings = loadings, scores = scores)
}
