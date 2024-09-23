#' @title Plot principal component confidence ellipses
#' @description Plot principal component confidence ellipses from multiple replicates
#'
#' @param results output of a \code{mambo} run.
#' @param locus label name of response or predictor locus.
#' @param pc.x number of x-axis principal component.
#' @param pc.y number of y-axis principal component.
#' @param type plot as ellipse of samples or 2-D density.
#' @param ellipse.p probability density level of ellipse.
#' @param num.bins number of bins for each axis if 2-D density is plotted.
#' @param plot display plot?
#'
#' @return PCA biplot of scores with confidence ellipses for each sample.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
plotPCs <- function(results, locus, pc.x = 1, pc.y = 2, 
                    type = c('ellipse', 'density'), ellipse.p = 0.95, 
                    num.bins = 50, plot = TRUE) {
  if(missing(locus)) stop("'locus' must be specified.")
  scores <- extractPCA(results)$scores[[locus]]
  type <- match.arg(type)
  
  df <- if(type == 'ellipse') {
    scores |> 
      split(scores$sample) |> 
      purrr::imap(function(df, i) {
        x <- dplyr::filter(df, pc == pc.x)$score
        y <- dplyr::filter(df, pc == pc.y)$score
        car::dataEllipse(x, y, levels = ellipse.p, draw = FALSE) |> 
          as.data.frame() |> 
          dplyr::mutate(sample = i)
      }) |> 
      dplyr::bind_rows()
  } else {
    scores |> 
      dplyr::mutate(axis = ifelse(pc == pc.x, 'x', 'y')) |> 
      dplyr::filter(pc %in% c(pc.x, pc.y)) |> 
      dplyr::select(-dplyr::all_of('pc')) |> 
      tidyr::pivot_wider(
        id_cols = c('sample', 'rep'), 
        names_from = 'axis', 
        values_from = 'score'
      )
  }
    
  gg <- df |> 
    ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_hline(yintercept = 0, color = 'darkred') +
    ggplot2::geom_vline(xintercept = 0, color = 'darkred') +
    ggplot2::labs(
      x = paste0('PC', pc.x),
      y = paste0('PC', pc.y),
      title = locus
    ) +
    ggplot2::theme_minimal()
  
  gg <- if(type == 'ellipse') {
    gg + 
      ggplot2::geom_polygon(
        ggplot2::aes(group = sample), 
        fill = NA, 
        color = 'black'
      ) 
  } else {
    gg +
      ggplot2::geom_bin_2d(bins = num.bins) +
      ggplot2::scale_fill_viridis_c(option = 'viridis') +
      ggplot2::theme(legend.position = 'none')
  }
  print(gg)
  
  if(plot) print(gg)
  invisible(gg)
}
