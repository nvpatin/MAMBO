#' @title Summary of principal component contrasts
#' @description Separate and sort positive and negative components
#'
#' @param results output of a \code{mambo} run.
#' @param d response or predictor label.
#' @param pc principal component to summarize (1:n).
#' @param min.rep minimum number of replicates to include.
#'
#' @return a list with the greatest positively (\code{$pos}) and
#'   negatively (\code{$neg}) loading ASVs sorted by magnitude.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
contrastSummary <- function(results, d, pc, min.rep = length(results)) {
  res <- lapply(results, function(x) {
    outlierLoadings(x$pca[[d]])[[pc]]
  })

  pos <- lapply(res, function(x) names(x$pos)) |>
    unlist() |>
    table() |>
    sort(decreasing = T)

  neg <- lapply(res, function(x) names(x$neg)) |>
    unlist() |>
    table() |>
    sort()

  list(
    pos = names(pos[pos >= min.rep]),
    neg = names(neg[neg >= min.rep])
  )
}
