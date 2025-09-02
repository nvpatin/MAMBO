#' @title Summarize PCA replicates
#' @description Number of components retained and distribution of cumulative 
#'   proportion of variance for by locus.
#'
#' @param results output of a \link[mambo]{mambo} run.
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#' 
#' @export
#'
pcaSummary <- function(results) {
  results$reps |> 
    lapply(function(x) x$pca) |> 
    purrr::list_transpose() |> 
    lapply(purrr::list_transpose) |>
    purrr::imap(function(reps, locus) {
      lapply(reps$importance, function(imp) {
        data.frame(
          locus = locus,
          num.pcs = ncol(imp),
          cum.prop = imp['Cumulative Proportion', ncol(imp)]
        )
      }) |> 
        dplyr::bind_rows() |> 
        dplyr::mutate(rep = 1:dplyr::n())
    }) |> 
    dplyr::bind_rows() |> 
    dplyr::select(dplyr::all_of(c('locus', 'rep', 'num.pcs', 'cum.prop')))
}
