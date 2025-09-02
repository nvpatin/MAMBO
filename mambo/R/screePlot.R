#' @title Scree plot of principal components
#' @description Plot median percent of variance accounted for by each 
#' principal component across multiple \code{MAMBO} replicates.
#'
#' @param results output of a \code{mambo} run.
#' @param plot display plot?
#'
#' @return \code{ggplot2} Scree plot of median percent of variance (y-axis) 
#' accounted for by each successive prinicpal component (x-axis). NB: Only
#' the minimum number of components occurring in all replicates are summarized.
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#'
#' @export
#'
screePlot <- function(results, plot = TRUE) {
  imp <- extractPCA(results)$importance
  p <- imp |> 
    dplyr::filter(.data$type %in% c('Proportion of Variance', 'Cumulative Proportion')) |> 
    dplyr::group_by(.data$locus, .data$type, .data$pc) |> 
    dplyr::summarize(
      value = stats::median(.data$value) * 100,
      .groups = 'drop'
    ) |> 
    ggplot2::ggplot(ggplot2::aes(x = .data$pc, y = .data$value)) +
    ggplot2::geom_line() +
    ggplot2::geom_point( 
      color = 'white', 
      fill = 'red', 
      shape = 21, 
      size = 4
    ) +
    ggplot2::scale_x_continuous(breaks = 1:max(imp$pc)) +
    ggplot2::labs(
      x = 'Principal Component', 
      y = 'Proportion'
    ) +
    ggplot2::facet_grid(.data$type ~ .data$locus, scales = 'free_y') +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(),
      axis.line.y = ggplot2::element_line()
    )
  
  if(plot) print(p)
  invisible(p)
}