#' @title Principal component analysis of random draw from beta distributon
#' @description Conducts a random draw of relative percent occurrence from 
#'   beta distribution. Then does a PCA of the log-odds of the relative percent
#'   occurence.
#'
#' @param beta.params matrix of beta parameters from \link{betaParams}.
#'
#' @return summary of PCA from \link{summary.prcomp}.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
ranPCA <- function(beta.params) {
  prob <- t(ranRelPct(beta.params))
  log(prob / (1 - prob)) |>
    stats::prcomp() |>
    summary()
}
