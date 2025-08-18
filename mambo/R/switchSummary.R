#' @title Summary of bayesian switches
#' @description Summarize posterior of predictor switches
#'
#' @param results output of a \link[mambo]{mambo} run.
#' @param min.p minimum proportion of inclusion to highlight a predictor.
#' @param resp.pcs,pred.pcs numeric vectors giving the response and predictor
#' PCs to summarize. If \code{NULL}, then the minimum number of PCs extracted 
#' for the response and predictor variables across all replicates are
#' summarized.
#' @param plot display summary plot?
#'
#' @return a list with a summary table and plot object.
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#'
#' @export
#'
switchSummary <- function(
    results, min.p = 0.9, resp.pcs = NULL, pred.pcs = NULL, plot = TRUE
) {
  resp <- results$labels['resp']
  pred <- results$labels['pred']

  min.pcs <- results$reps |>
    sapply(function(x) {
      sapply(x$pca, function(pca.x) ncol(pca.x$x))
    }) |>
    apply(1, min)
  
  resp.pcs <- if(is.null(resp.pcs)) 1:min.pcs[resp] else {
    pcs <- 1:min.pcs[resp]
    pcs <- resp.pcs[resp.pcs %in% pcs]
    if(length(pcs) == 0) {
      stop('Requested response PCs not available.')
    }
    pcs
  }
    
  pred.pcs <- if(is.null(pred.pcs)) 1:min.pcs[pred] else {
    pcs <- 1:min.pcs[pred]
    pcs <- pred.pcs[pred.pcs %in% pcs]
    if(length(pcs) == 0) {
      stop('Requested predictor PCs not available.')
    }
    pcs
  }
  
  w.post <- do.call(
    abind::abind,
    c(lapply(results$reps, function(x) {
      apply(x$post.list$w[resp.pcs, pred.pcs, ], c(1, 2), mean)
    }), list(along = 3))
  )
  dimnames(w.post)[[3]] <- 1:dim(w.post)[3]
  names(dimnames(w.post)) <- c(resp, pred, 'rep')

  w.post <- w.post |>
    as.data.frame.table() |>
    stats::setNames(c('resp', 'pred', 'rep', 'w')) |>
    dplyr::mutate(rep = as.numeric(.data$rep))

  smry <- w.post |>
    dplyr::group_by(resp, pred) |>
    dplyr::summarize(median.w = stats::median(.data$w), .groups = 'drop') |>
    dplyr::filter(.data$median.w > min.p)

  p <- w.post |>
    dplyr::left_join(smry, by = c('resp', 'pred')) |>
    dplyr::mutate(to.highlight = !is.na(.data$median.w)) |>
    ggplot2::ggplot() +
    ggplot2::geom_histogram(
      ggplot2::aes(.data$w, fill = .data$to.highlight),
      binwidth = 0.05
    ) +
    ggplot2::xlim(c(0, 1)) +
    ggplot2::scale_fill_manual(values = c('black', 'red')) +
    ggplot2::facet_grid(pred ~ resp) +
    ggplot2::theme(
      legend.position = 'none',
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  if(plot) print(p)

  list(smry = smry, plot = p)
}
