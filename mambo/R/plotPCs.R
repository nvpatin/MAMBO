#' @title Plot principal component confidence ellipses
#' @description Plot principal component confidence ellipses from multiple replicates
#'
#' @param results output of a \code{mambo} run.
#' @param locus label name of response or predictor locus.
#' @param pc.x number of x-axis principal component.
#' @param pc.y number of y-axis principal component.
#' @param ellipse.p probability density level of ellipse.
#'
#' @return PCA biplot of scores with confidence ellipses for each sample.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
plotPCs <- function(results, locus, pc.x = 1, pc.y = 2, ellipse.p = 0.95) {
  scores <- extractPCA(results)$scores[[locus]]
  
  purrr::imap(split(scores, scores$sample), function(df, i) {
    x <- dplyr::filter(df, pc == pc.x)$score
    y <- dplyr::filter(df, pc == pc.y)$score
    car::dataEllipse(x, y, levels = ellipse.p, draw = FALSE) |> 
      as.data.frame() |> 
      dplyr::mutate(sample = i)
  }) |> 
    dplyr::bind_rows() |> 
    ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, color = 'darkred') +
    ggplot2::geom_vline(xintercept = 0, color = 'darkred') +
    ggplot2::geom_polygon(
      ggplot2::aes(x, y, group = sample),
      fill = NA,
      color = 'black'
    ) +
    ggplot2::labs(
      x = paste0('PC', pc.x),
      y = paste0('PC', pc.y)
    ) +
    theme_minimal()
}
