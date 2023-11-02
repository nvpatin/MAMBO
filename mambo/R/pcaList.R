#' @title Extract transposed list of all PCA results.
#' @description Extract transposed list of all PCA results.
#'
#' @param results output of a \code{mambo} run.
#'
#' @return a list of PCA results by type across all replicates.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
pcaList <- function(results) {
  results$reps |> 
    lapply(function(x) x$pca) |> 
    purrr::list_transpose() |> 
    lapply(purrr::list_transpose)
}