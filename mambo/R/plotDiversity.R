#' @title Plot diversity of predictors and response loci
#' @description Plot diversity (effective number of ASVs) for predictor and
#' response loci.
#'
#' @param results output of a \code{mambo} run.
#' @param ellipse.p probability density level of ellipse.
#' @param sample.df a data frame containing columns to color ellipses by 
#' (\code{ellipse.fill}) or facet by (\code{facet.by}). Must have a column 
#' called '\code{sample}' identifying samples in \code{locus}.
#' @param ellipse.fill a column in the \code{sample.df} data frame to use to 
#' color in ellipses.
#' @param facet.by a column in the \code{sample.df} data frame to use to facet 
#' the plot by.
#' @param log.axes log-scale axes?
#' @param plot display plot?
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#'
#' @export
#'
plotDiversity <- function(
    results, ellipse.p = 0.95,  sample.df = NULL, ellipse.fill = NULL, 
    facet.by = NULL, log.axes = TRUE,  plot = TRUE
) {
  if(is.na(results$labels['pred'])) {
    stop("'results' must have both response and predictor values.")
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
  
  resp <- results$labels['resp']
  pred <- results$labels['pred']
  
  df <- results$reps |> 
    lapply(function(x) x$sample.diversity) |> 
    dplyr::bind_rows() 
  
  df <- split(df, df$sample) |> 
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
  
  if(!is.null(sample.df)) df <- dplyr::left_join(df, sample.df, by = 'sample')
  
  gg <- df |> 
    ggplot2::ggplot(mapping = ggplot2::aes(x = .data[[pred]], y = .data[[resp]])) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = 'darkred') +
    ggplot2::labs(
      x = paste0('Effective Number of ASVs (', pred, ')'),
      y = paste0('Effective Number of ASVs (', resp, ')')
    ) +
    ggplot2::theme_minimal()
  
  gg <- if(!is.null(ellipse.fill)) {
    x <- gg + ggplot2::geom_polygon(
      ggplot2::aes(group = .data$sample, fill = .data[[ellipse.fill]]),
      alpha = 0.8,
      color = 'white'
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
  
  if(log.axes) {
    gg <- gg +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10()
  }
  
  if(!is.null(facet.by)) {
    gg <- gg + ggplot2::facet_wrap(~ .data[[facet.by]], ncol = 1)
  }
  
  if(plot) print(gg)
  invisible(gg)
}
