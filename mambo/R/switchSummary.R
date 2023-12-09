#' @title Summary of bayesian switches
#' @description Summarize posterior of predictor switches
#'
#' @param results output of a \code{mambo} run.
#' @param min.p minimum proportion of inclusion to highlight a predictor.
#' @param plot display summary plot?
#'
#' @return a list with a summary table and plot object.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
switchSummary <- function(results, min.p = 0.75, plot = TRUE) {
  resp <- results$labels['resp']
  pred <- results$labels['pred']

  min.pcs <- results$reps |>
    sapply(function(x) {
      sapply(x$pca, function(pca.x) pca.x$num.pcs)
    }) |>
    apply(1, min)

  w.post <- do.call(
    abind::abind,
    c(lapply(results$reps, function(x) {
      apply(x$post.list$w[1:min.pcs[resp], 1:min.pcs[pred], ], c(1, 2), mean)
    }), list(along = 3))
  )
  dimnames(w.post)[[3]] <- 1:dim(w.post)[3]
  names(dimnames(w.post)) <- c(resp, pred, 'rep')

  w.post <- w.post |>
    as.data.frame.table() |>
    stats::setNames(c('resp', 'pred', 'rep', 'w')) |>
    dplyr::mutate(rep = as.numeric(rep))

  smry <- w.post |>
    dplyr::group_by(resp, pred) |>
    dplyr::summarize(median.w = stats::median(w), .groups = 'drop') |>
    dplyr::filter(median.w > min.p)

  p <- w.post |>
    dplyr::left_join(smry, by = c('resp', 'pred')) |>
    dplyr::mutate(to.highlight = !is.na(median.w)) |>
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(w, fill = to.highlight)) +
    ggplot2::scale_fill_manual(values = c('black', 'red')) +
    ggplot2::facet_grid(pred ~ resp) +
    ggplot2::theme(legend.position = 'none')

  if(plot) print(p)

  list(smry = smry, plot = p)
}
