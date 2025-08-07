#' @title Plot diversity of predictors and response loci
#' @description Plot diversity (effective number of ASVs) for predictor and
#' response loci
#'
#' @param results output of a \code{mambo} run.
#' @param type plot as ellipse of samples or 2-D density.
#' @param ellipse.p probability density level of ellipse.
#' @param num.bins number of bins for each axis if 2-D density is plotted.
#' @param sample.df a data frame containing columns to color ellipses by 
#' (\code{ellipse.fill}) or facet by (\code{facet.by}). Must have a column 
#' called '\code{sample}' identifying samples in \code{locus}.
#' @param ellipse.fill a column in the \code{sample.df} data frame to use to 
#' color in ellipses.
#' @param facet.by a column in the \code{sample.df} data frame to use to facet 
#' the plot by.
#' @param plot display plot?
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
plotDiversity <- function(results, type = c('ellipse', 'density'), 
                          ellipse.p = 0.95, num.bins = 50, sample.df = NULL, 
                          ellipse.fill = NULL, facet.by = NULL, plot = TRUE) {
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

  df <- results$reps |> 
    lapply(function(x) x$diversity) |> 
    dplyr::bind_rows()
  
  pred <- results$labels['pred']
  resp <- results$labels['resp']
  
  type <- match.arg(type)
  if(type == 'ellipse') {
    df <- df |> 
      split(df$sample) |> 
      purrr::imap(function(df, i) {
        car::dataEllipse(
          df[[pred]], 
          df[[resp]], 
          levels = ellipse.p, 
          draw = FALSE
        ) |> 
          as.data.frame() |> 
          dplyr::mutate(sample = i)
      }) |> 
      dplyr::bind_rows() |> 
      stats::setNames(c(pred, resp, 'sample'))
  }
  
  if(!is.null(sample.df)) df <- dplyr::left_join(df, sample.df, by = 'sample')
    
  gg <- df |> 
    ggplot2::ggplot(mapping = ggplot2::aes(x = .data[[pred]], y = .data[[resp]])) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = 'darkred') +
    ggplot2::labs(
      x = paste0('Effective Number of ASVs (', pred, ')'),
      y = paste0('Effective Number of ASVs (', resp, ')')
    ) +
    ggplot2::theme_minimal()
  
  gg <- if(type == 'ellipse') {
    if(!is.null(ellipse.fill)) {
      gg + ggplot2::geom_polygon(
        ggplot2::aes(group = .data$sample, fill = .data[[ellipse.fill]]),
        alpha = 0.75,
        color = 'black'
      ) +
        ggplot2::scale_fill_viridis_c(option = 'viridis') +
        ggplot2::theme(legend.position = 'top')
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
