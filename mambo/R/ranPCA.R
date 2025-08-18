#' @title Principal component analysis of random draw from beta distribution
#' @description Computes a PCA of the log-odds of the relative percent
#'   occurrence.
#'
#' @param rel.abund matrix of relative abundance expressed as percents of 
#' read counts (\code{0:1}). Columns are samples, rows are ASVs.
#'
#' @return Summary of PCA from \link{summary.prcomp}.
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#'
ranPCA <- function(rel.abund) {
  prob <- t(rel.abund)
  log(prob / (1 - prob)) |>
    stats::prcomp() |>
    summary() |> 
    suppressWarnings()
}
