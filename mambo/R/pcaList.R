#' @title Extract transposed list of all PCA results.
#' @description Extract transposed list of all PCA results.
#'
#' @param results output of a \link[mambo]{mambo} run.
#'
#' @return a list of PCA results by type across all replicates.
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#'
#' @export
#'
pcaList <- function(results) {
  results$reps |> 
    lapply(function(x) x$pca) |> 
    purrr::list_transpose() |> 
    lapply(purrr::list_transpose)
}
