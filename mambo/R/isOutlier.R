#' @title Identify outlier values
#' @description Separate and sort positive and negative components
#'
#' @param x a vector of numbers.
#' @param type choose outlier based on inter-quartile interval (\code{iqr}) or
#'   z-score (\code{z}).
#' @param thresh if \code{type = 'z'}, the minimum z-score used to identify
#'   outliers.
#'
#' @return a vector of logicals identifying values that are outliers.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
isOutlier <- function(x, type = c('iqr', 'z'), thresh = 3) {
  switch(
    match.arg(type),
    iqr = {
      quarts <- stats::quantile(x, probs = c(0.25, 0.75))
      iqr <- diff(quarts)
      thresh <- c(quarts[1] - 1.5 * iqr, quarts[2] + 1.5 * iqr)
      x <= thresh[1] | x >= thresh[2]
    },
    z = {
      z.score <- (x - mean(x)) / stats::sd(x)
      abs(z.score) >= thresh
    },
    NULL
  )
}
