#' @title Extract PCA values
#' @description Extract PCA loadings, scores, and importance metrics.
#'
#' @param results output of a \code{mambo} run.
#'
#' @return a list with data frames of component loadings (\code{$loadings}),
#' scores (\code{$scores}), and component importance metrics (\code{$importance})
#' for each replicate and locus. 
#' NB: Only values for the minimum number of components across all replicates
#' are returned.
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#'
#' @export
#'
extractPCA <- function(results) {
  pca.list <- results$reps |> 
    lapply(function(x) x$pca) |> 
    purrr::list_transpose() |> 
    lapply(purrr::list_transpose)
  
  .switchSigns <- function(x) {
    # use sample (row) with largest absolute loading on first replicate (column) as reference for sign
    ref.index <- which.max(abs(x[, 1]))
    # identify replicates (columns) that have a different sign from reference
    to.switch <- sign(x[ref.index, ]) != sign(x[ref.index, 1])
    # switch sign of replicates that need it
    x[, to.switch] <- x[, to.switch] * -1
    x
  }
  
  loadings <- purrr::imap(pca.list, function(reps, locus) {
    min.pcs <- min(sapply(reps$rotation, function(x) dim(x)[2]))
    res <- reps$rotation |> 
      lapply(function(x) x[, 1:min.pcs]) |> 
      abind::abind(along = 3) |> 
      apply(2, .switchSigns, simplify = FALSE) |> 
      abind::abind(along = 3) |> 
      aperm(c(1, 3, 2))
    dimnames(res)[[3]] <- 1:dim(res)[3]
    res |> 
      as.data.frame.table() |> 
      stats::setNames(c('asv', 'pc', 'rep', 'loading')) |> 
      dplyr::mutate(
        locus = locus,
        asv = as.character(.data$asv),
        pc = as.numeric(stringr::str_remove(.data$pc, 'PC')),
        rep = as.numeric(.data$rep)
      )
  }) |> 
    dplyr::bind_rows()
  
  scores <- purrr::imap(pca.list, function(reps, locus) {
    min.pcs <- min(sapply(reps$x, function(x) dim(x)[2]))
    res <- reps$x |> 
      lapply(function(x) x[, 1:min.pcs]) |> 
      abind::abind(along = 3) |> 
      apply(2, .switchSigns, simplify = FALSE) |> 
      abind::abind(along = 3) |> 
      aperm(c(1, 3, 2))
    dimnames(res)[[3]] <- 1:dim(res)[3]
    res |> 
      as.data.frame.table() |> 
      stats::setNames(c('sample', 'pc', 'rep', 'score')) |> 
      dplyr::mutate(
        locus = locus,
        sample = as.character(.data$sample),
        pc = as.numeric(stringr::str_remove(.data$pc, 'PC')),
        rep = as.numeric(.data$rep)
      )
  }) |> 
    dplyr::bind_rows()
  
  importance <- purrr::imap(pca.list, function(reps, locus) {
    min.pcs <- min(sapply(reps$importance, function(x) dim(x)[2]))
    res <- reps$importance |> 
      lapply(function(x) x[, 1:min.pcs]) |> 
      abind::abind(along = 3)
    dimnames(res)[[3]] <- 1:dim(res)[3]
    res <- res |> 
      as.data.frame.table() |> 
      stats::setNames(c('type', 'pc', 'rep', 'value')) |> 
      dplyr::mutate(
        locus = locus,
        type = as.character(.data$type),
        pc = as.numeric(stringr::str_remove(.data$pc, 'PC')),
        rep = as.numeric(.data$rep)
      )
  }) |> 
    dplyr::bind_rows()
  
  list(loadings = loadings, scores = scores, importance = importance)
}
