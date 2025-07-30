#' @title Summarize PCA replicates
#' @description Number of components retained and distribution of cumulative 
#'   proportion of variance for by locus.
#'
#' @param results output of a \link[mambo]{mambo} run.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#'
pcaSummary <- function(results) {
  lapply(results$reps, function(x) {
    sapply(x$pca, function(pca) {
      num.pcs <- ncol(pca$importance)
      c(num.pcs = num.pcs, cum.pct.var = pca$importance[3, num.pcs])
    }, USE.NAMES = TRUE) |> 
      t() |> 
      as.data.frame() |> 
      tibble::rownames_to_column('locus') 
  }) |> 
    dplyr::bind_rows()
}
