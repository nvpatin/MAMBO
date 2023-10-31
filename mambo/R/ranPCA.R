#' @title Principal component analysis of random draw from beta distributon
#' @description Principal component analysis of random draw from beta distributon
#'
#' @param beta.params matrix of beta parameters from \link{betaParams}.
#'
#' @return summary list of PCA of log-odds of random draw.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
ranPCA <- function(beta.params) {
  prob <- t(ranRelPct(beta.params))
  log(prob / (1 - prob)) |>
    stats::prcomp() |>
    summary()
}
