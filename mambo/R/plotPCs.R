#' @title Plot principal component scores with replicate uncertainty
#' @description Plot principal component scores as ellipses or 2-D densities
#' from multiple \code{MAMBO} replicates.
#'
#' @param results output of a \code{mambo} run.
#' @param locus label name of response or predictor locus.
#' @param pc.x number of x-axis principal component.
#' @param pc.y number of y-axis principal component.
#' @param type plot as ellipse of samples or 2-D density.
#' @param ellipse.p probability density level of ellipse.
#' @param num.bins number of bins for each axis if 2-D density is plotted.
#' @param sample.df a data frame containing columns to color ellipses by 
#' (\code{ellipse.fill}) or facet by (\code{facet.by}). Must have a column 
#' called '\code{sample}' identifying samples in \code{locus}.
#' @param ellipse.fill a column in the \code{sample.df} data frame to use to 
#' color in ellipses. If this column is a \code{numeric} it is treated as a
#' continuous variable, otherwise it is treated as discrete.
#' @param facet.by a column in the \code{sample.df} data frame to use to facet 
#' the plot by.
#' @param plot display plot?
#'
#' @return \code{ggplot2} PCA biplot of scores.
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#'
#' @export
#'
plotPCs <- function(results, locus, pc.x = 1, pc.y = 2, 
                    type = c('ellipse', 'density'), ellipse.p = 0.95, 
                    num.bins = 50, sample.df = NULL, ellipse.fill = NULL,
                    facet.by = NULL, plot = TRUE) {
  if(missing(locus)) stop("'locus' must be specified.")
  if(!locus %in% results$labels[c('resp', 'pred')]) {
    stop("locus '", locus, "' is not in 'results'")
  }
  
  if(!is.null(ellipse.fill)) {
    if(!ellipse.fill %in% names(sample.df)) {
      stop("'", ellipse.fill, "' not in 'sample.df'")
    }
  }
  
  if(!is.null(facet.by)) {
    if(!facet.by %in% names(sample.df)) {
      stop("'", facet.by, "' not in 'sample.df'")
    }
  }
  
  scores <- extractPCA(results)$scores[[locus]]
  
  prop.var <- sapply(
    results$reps, 
    function(r) r$pca[[locus]]$importance['Proportion of Variance', c(pc.x, pc.y)]
  ) |> 
    t() |> 
    apply(2, stats::median)

  type <- match.arg(type)
  df <- if(type == 'ellipse') {
    scores |> 
      split(scores$sample) |> 
      purrr::imap(function(df, i) {
        x <- dplyr::filter(df, .data$pc == pc.x)$score
        y <- dplyr::filter(df, .data$pc == pc.y)$score
        car::dataEllipse(x, y, levels = ellipse.p, draw = FALSE) |> 
          as.data.frame() |> 
          dplyr::mutate(sample = i)
      }) |> 
      dplyr::bind_rows()
  } else {
    scores |> 
      dplyr::mutate(axis = ifelse(.data$pc == pc.x, 'x', 'y')) |> 
      dplyr::filter(.data$pc %in% c(pc.x, pc.y)) |> 
      dplyr::select(-dplyr::all_of('pc')) |> 
      tidyr::pivot_wider(
        id_cols = c('sample', 'rep'), 
        names_from = 'axis', 
        values_from = 'score'
      )
  }
  
  if(!is.null(sample.df)) df <- dplyr::left_join(df, sample.df, by = 'sample')
    
  gg <- df |> 
    ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_hline(yintercept = 0, color = 'darkred') +
    ggplot2::geom_vline(xintercept = 0, color = 'darkred') +
    ggplot2::labs(
      x = paste0('PC', pc.x, ' (', round(prop.var[1] * 100, 1), '%)'),
      y = paste0('PC', pc.y, ' (', round(prop.var[2] * 100, 1), '%)'),
      title = locus
    ) +
    ggplot2::theme_minimal()
  
  gg <- if(type == 'ellipse') {
    if(!is.null(ellipse.fill)) {
      x <- gg + ggplot2::geom_polygon(
        ggplot2::aes(group = .data$sample, fill = .data[[ellipse.fill]]),
        alpha = 0.75,
        color = 'black'
      ) 
      if(is.numeric(df[[ellipse.fill]])) {
        x <- x + ggplot2::scale_fill_viridis_c(option = 'viridis')
      } else {
        x <- x + ggplot2::scale_fill_brewer(palette = 'Set2')
      }
      x + ggplot2::theme(legend.position = 'top')
    } else {
      gg + ggplot2::geom_polygon(
        ggplot2::aes(group = .data$sample), 
        fill = NA, 
        color = 'black'
      )
    }
  } else {
    gg +
      ggplot2::geom_bin_2d(bins = num.bins) +
      ggplot2::scale_fill_viridis_c(option = 'viridis') +
      ggplot2::theme(legend.position = 'none')
  }
    
  if(!is.null(facet.by)) gg <- gg + ggplot2::facet_wrap(~ .data[[facet.by]], ncol = 1)
  
  if(plot) print(gg)
  invisible(gg)
}
