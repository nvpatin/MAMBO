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
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
outlierLoadings <- function(
  results, locus, pc, type = c('iqr', 'z'), thresh = 3, min.pct.reps = 0.95
) {
  extractPCA(results)$loadings[[locus]] |> 
    dplyr::rename(pcs = 'pc') |> 
    dplyr::filter(pcs == pc) |> 
    dplyr::group_by(rep) |> 
    dplyr::mutate(is.outlier = isOutlier(loading, type, thresh)) |> 
    dplyr::ungroup() |> 
    dplyr::group_by(asv) |> 
    dplyr::summarize(
      mean.loading = mean(loading),
      median.loading = stats::median(loading),
      pct.reps = mean(is.outlier),
      .groups = 'drop'
    ) |> 
    dplyr::arrange(dplyr::desc(median.loading)) |> 
    dplyr::filter(pct.reps >= min.pct.reps)
}
