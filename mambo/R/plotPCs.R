#' @title Plot principal component scores with replicate uncertainty
#' @description Create biplot of principal component scores as ellipses
#' from multiple \code{MAMBO} replicates.
#'
#' @param results output of a \code{mambo} run.
#' @param locus.x label name of response or predictor locus for x-axis. If not 
#' specified, then it defaults to the response locus. 
#' @param locus.y label name of response or predictor locus fo y-axis. If not 
#' specified, then it defaults to \code{locus.x}. The plot will then be a 
#' PCA biplot for a single locus.
#' @param pc.x number of principal component for x-axis.
#' @param pc.y number of principal component for y-axis.
#' @param ellipse.p probability density level of ellipse.
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
plotPCs <- function(
    results, locus.x = NULL, locus.y = locus.x, pc.x = 1, pc.y = 2, 
    ellipse.p = 0.95, sample.df = NULL, ellipse.fill = NULL,
    facet.by = NULL, plot = TRUE
) {
  if(is.null(locus.x)) locus.x <- results$labels['resp']
  
  if(!locus.x %in% results$labels[c('resp', 'pred')]) {
    stop("locus '", locus.x, "' is not in 'results'")
  }
  if(!locus.y %in% results$labels[c('resp', 'pred')]) {
    stop("locus '", locus.y, "' is not in 'results'")
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
  
  pc.list <- extractPCA(results)
  
  prop.var <- pc.list$importance |> 
    dplyr::filter(
      ((.data$locus == locus.x & .data$pc == pc.x) |
         (.data$locus == locus.y & .data$pc == pc.y)) &
        .data$type == 'Proportion of Variance'
    ) |> 
    dplyr::group_by(.data$locus, .data$pc) |> 
    dplyr::summarize(
      prop = round(stats::median(.data$value) * 100, 1),
      .groups = 'drop'
    ) |> 
    dplyr::mutate(
      axis = ifelse(.data$locus == locus.x & .data$pc == pc.x, 'x', 'y'),
      label = paste0(.data$locus, ' PC', .data$pc, ' (', .data$prop, '%)')
    ) |> 
    dplyr::arrange(.data$axis) |> 
    dplyr::select(dplyr::all_of(c('axis', 'label'))) |> 
    tibble::deframe()
  
  df <- pc.list$scores |> 
    split(pc.list$scores$sample) |> 
    purrr::imap(function(df, i) {
      x <- dplyr::filter(df, .data$pc == pc.x & .data$locus == locus.x)$score
      y <- dplyr::filter(df, .data$pc == pc.y & .data$locus == locus.y)$score
      car::dataEllipse(x, y, levels = ellipse.p, draw = FALSE) |> 
        as.data.frame() |> 
        dplyr::mutate(sample = i)
    }) |> 
    dplyr::bind_rows() 
  
  if(!is.null(sample.df)) df <- dplyr::left_join(df, sample.df, by = 'sample')
  
  gg <- df |> 
    ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_hline(yintercept = 0, color = 'darkred') +
    ggplot2::geom_vline(xintercept = 0, color = 'darkred') +
    ggplot2::labs(x = prop.var['x'], y = prop.var['y']) +
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
  
  if(!is.null(facet.by)) gg <- gg + ggplot2::facet_wrap(~ .data[[facet.by]], ncol = 1)
  
  if(plot) print(gg)
  invisible(gg)
}