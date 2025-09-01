#' @title Scree plot of principal components
#' @description Plot median percent of variance accounted for by each 
#' principal component across multiple \code{MAMBO} replicates.
#'
#' @param results output of a \code{mambo} run.
#' @param locus label name of response or predictor locus.
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
screePlot <- function(results, locus, plot = TRUE) {
  min.pc <- min(
    sapply(results$reps, function(r) ncol(r$pca[[locus]]$importance))
  )
  
  p <- sapply(
    results$reps, 
    function(r) r$pca[[locus]]$importance['Proportion of Variance', 1:min.pc]
  ) |> 
    t() |> 
    as.data.frame() |> 
    dplyr::mutate(rep = 1:dplyr::n()) |> 
    tidyr::pivot_longer(-rep, names_to = 'pc', values_to = 'pct') |> 
    dplyr::mutate(pct = .data$pct * 100) |> 
    dplyr::group_by(.data$pc) |> 
    dplyr::summarize(pct.var = stats::median(.data$pct), .groups = 'drop') |> 
    dplyr::mutate(pc = as.numeric(stringr::str_remove(.data$pc, 'PC'))) |> 
    dplyr::arrange(.data$pc) |> 
    dplyr::mutate(cum.pct = cumsum(pct.var)) |> 
    tidyr::pivot_longer(-pc, names_to = 'type', values_to = 'pct') |> 
    dplyr::mutate(
      type = ifelse(type == 'pct.var', 'Absolute', 'Cumulative')
    ) |> 
    ggplot2::ggplot(ggplot2::aes(x = .data$pc)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$pct)) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$pct), 
      color = 'white', 
      fill = 'red', 
      shape = 21, 
      size = 4
    ) +
    ggplot2::scale_x_continuous(breaks = 1:min.pc) +
    ggplot2::labs(
      x = 'Principal Component', 
      y = 'Percent of Variance',
      title = locus
    ) +
    ggplot2::facet_wrap(~ type, nrow = 2, scales = 'free_y') +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(),
      axis.line.y = ggplot2::element_line()
    )
  
  if(plot) print(p)
  invisible(p)
}