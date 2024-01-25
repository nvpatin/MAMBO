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
  pca.list <- pcaList(results)
  
  loadings <- lapply(pca.list, function(rep) {
    res <- abind::abind(rep$rotation, along = 3) |> 
      apply(2, function(x) {
        # use ASV (row) with largest absolute loading on first replicate (column) as reference for sign
        ref.index <- which.max(abs(x[, 1]))
        # identify replicates (columns) that have a different sign from reference
        to.switch <- sign(x[ref.index, ]) != sign(x[ref.index, 1])
        # switch sign of replicates that need it
        x[, to.switch] <- x[, to.switch] * -1
        x
      }, simplify = FALSE) |> 
      abind::abind(along = 3) |> 
      aperm(c(1, 3, 2))
    dimnames(res)[[3]] <- 1:dim(res)[3]
    as.data.frame.table(res) |> 
      stats::setNames(c('asv', 'pc', 'rep', 'loading')) |> 
      dplyr::mutate(
        asv = as.character(asv),
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
