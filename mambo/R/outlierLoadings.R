#' @title Outlier principal component loadings
#' @description Return loading summary of outlier ASVs.
#'
#' @param results output of a \link{mambo} run.
#' @param locus label name of response or predictor locus.
#' @param pc principal component to summarize (1:n).#' 
#' @param type choose outlier based on inter-quartile interval (\code{iqr}) or
#'   z-score (\code{z}).
#' @param thresh if \code{type = 'z'}, the minimum z-score used to identify
#'   outliers.
#' @param min.pct.reps minimum percent of replicates that an ASV should be 
#'   identified as an outlier in to be saved.
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#'
#' @export
#'
outlierLoadings <- function(
  results, locus = NULL, pc = 1, type = c('z', 'iqr'), thresh = 3, min.pct.reps = 0.95
) {
  if(is.null(locus)) locus <- results$labels['resp']
  
  extractPCA(results)$loadings |> 
    dplyr::rename(p = 'pc', l = 'locus') |> 
    dplyr::filter(.data$p == pc & .data$l == locus) |> 
    dplyr::group_by(.data$rep) |> 
    dplyr::mutate(is.outlier = isOutlier(.data$loading, type, thresh)) |> 
    dplyr::ungroup() |> 
    dplyr::group_by(.data$asv) |> 
    dplyr::summarize(
      mean.loading = mean(.data$loading),
      median.loading = stats::median(.data$loading),
      pct.reps = mean(.data$is.outlier),
      .groups = 'drop'
    ) |> 
    dplyr::arrange(dplyr::desc(.data$median.loading)) |> 
    dplyr::filter(.data$pct.reps >= min.pct.reps)
}
