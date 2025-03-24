#' @title Principal component analysis of random draw from beta distribution
#' @description Conducts a random draw of relative percent occurrence from 
#'   beta distribution. Then does a PCA of the log-odds of the relative percent
#'   occurrence.
#'
#' @param beta.params matrix of beta parameters from \link{betaParams}.
#' @param num.cores number of cores to use for sampling relative percent 
#'    occurrence. Defaults to value reported by \link[parallel]{detectCores} 
#'    minus 1.
#'
#' @return summary of PCA from \link{summary.prcomp}.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
ranPCA <- function(beta.params, num.cores = parallel::detectCores() - 1) {
  prob <- t(ranRelPct(beta.params, num.cores))
  log(prob / (1 - prob)) |>
    stats::prcomp() |>
    summary() |> 
    suppressWarnings()
}