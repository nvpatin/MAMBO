#' @title Run MAMBO replicates
#' @description Run MAMBO replicates
#'
#' @param resp.label label for response locus.
#' @param resp.counts ASV counts of response locus.
#' @param pred.label label for predictor locus.
#' @param pred.counts ASV counts of predictor locus.
#' @param nrep number of MAMBO replicates to run.
#' @param chains number of MCMC chains.
#' @param adapt number of adaptation iterations.
#' @param burnin number of burnin iterations.
#' @param total.samples total number of samples from the posterior to save.
#' @param thin number of iterations to skip between samples in each chain.
#' @param run.label label for the run output.
#' @param output.log create a text log of the run progress?
#'
#' @note \code{resp.counts} and \code{pred.counts} should be matrices, 
#' data frames, or names of comma-delimited (.csv) files where values are the 
#' occurence (number of reads) for each ASV (rows) in each sample (columns).
#' 
#' @return a list containing:
#' \describe{ 
#'   \item{\code{$filename}}{the name of the RDS file written.}
#'   \item{\code{$labels}}{the run, response and predictor labels.} 
#'   \item{\code{$params}}{the MCMC run parameters.}
#'   \item{\code{$reps}}{a list with results for each of \code{nrep} replicate runs.}
#'   \item{\code{$run.time}}{a list of the start, end, and elapsed run times for all replicates.}
#' }
#' Each replicate in the \code{$reps} element contains:
#' \describe{
#'   \item{\code{$pca}}{a list of the response and predictor PCA results.}
#'   \item{\code{$post.smry}}{summary of parameter posterior for the replicate.}
#'   \item{\code{$post.list}}{list of posterior distribution for each parameter in the model.}
#' }
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
mambo <- function(
    resp.label, resp.counts, 
    pred.label, pred.counts,
    nrep = 10,
    chains = 3,
    adapt = 100,
    burnin = 1000,
    total.samples = 1000,
    thin = 1,
    run.label = 'mambo',
    output.log = TRUE
) {
  start.time <- Sys.time()
  if(output.log) {
    log.fname <- paste0(run.label, '.', format(start.time, '%Y%m%d_%H%M%S.log'))
    log.file <- file(log.fname, open = 'a')
    sink(file = log.file, type = 'output', split = TRUE)
  }
  
  cat('\n--------', format(Sys.time()), 'Starting MAMBO --------\n')
  cat('  Number of replicates:', nrep, '\n')
  cat('  MCMC parameters:\n')
  cat('    Chains:', chains, '\n')
  cat('    Adapt:', adapt, '\n')
  cat('    Burnin:', burnin, '\n')
  cat('    Total Samples:', total.samples, '\n')
  cat('    Thinning:', thin, '\n')
  if(output.log) cat('  Log file:', log.fname, '\n')
  
  cat('\n--------', format(Sys.time()), 'Occurrence Beta parameters --------\n')
  # compute beta parameter arrays
  resp.beta <- betaParams(resp.counts)
  pred.beta <- betaParams(pred.counts)
  
  # check that the same sample names are in both datasets
  if(!setequal(
    sort(dimnames(resp.beta)[[1]]), 
    sort(dimnames(pred.beta)[[1]])
  )) stop("sample names in 'resp.counts' and 'pred.counts' are not the same.")
  
  # make sure rows are in same order for both sets of data
  resp.beta <- resp.beta[dimnames(pred.beta)[[1]], , ]
  
  # do nrep iterations, save results to list, and write to RDS file
  reps <- lapply(1:nrep, function(i) {
    cat('\n--------', format(Sys.time()), 'Replicate ')
    cat(i, '/', nrep, sep = '')
    cat(' --------\n')
    
    # Extract PCs -------------------------------------------------------------
    cat('  PCA...\n')
    pca <- stats::setNames(
      list(ranPCA(resp.beta), ranPCA(pred.beta)),
      c(resp.label, pred.label)
    )
    pca[[resp.label]]$num.pcs <- numImpPCs(pca[[resp.label]])
    pca[[pred.label]]$num.pcs <- numImpPCs(pca[[pred.label]])
    
    if(bayesian) {
      # Run Bayesian model ------------------------------------------------------
      cat('  Bayesian model...\n')
      capture.output(post <- jagsPClm(
        pc.resp = pca[[resp.label]]$x[, 1:pca[[resp.label]]$num.pcs],
        pc.preds = pca[[pred.label]]$x[, 1:pca[[pred.label]]$num.pcs],
        chains = chains,
        adapt = adapt,
        burnin = burnin,
        total.samples = total.samples,
        thin = thin
      ))
      
      # Compute posterior summary statistics ------------------------------------
      cat('  Summarize posterior...\n')
      capture.output(post.smry <- summary(post, silent.jags = TRUE))
      
      # Extract posterior and label dimensions -----------------------------------
      p <- swfscMisc::runjags2list(post)
      dimnames(p$intercept)[[1]] <-
        dimnames(p$b.prime)[[1]] <-
        dimnames(p$w)[[1]] <-
        dimnames(p$v)[[1]] <- paste0(resp.label, '.PC', 1:pca[[resp.label]]$num.pcs)
      dimnames(p$b.prime)[[2]] <-
        dimnames(p$w)[[2]] <- paste0(pred.label, '.PC', 1:pca[[pred.label]]$num.pcs)
    }
    
    cat(
      '  End replicate:', 
      format(round(swfscMisc::autoUnits(post$timetaken))),
      '\n'
    )
    list(pca = pca, post.smry = post.smry, post.list = p)
  })
  end.time = Sys.time()
  
  res <- list(
    filename = paste0(run.label, '.', format(start.time, '%Y%m%d_%H%M%S.rds')),
    labels = c(run = run.label, resp = resp.label, pred = pred.label),
    params = c(
      nrep = nrep, chains = chains, adapt = adapt,
      burnin = burnin, total.samples = total.samples, thin = thin
    ),
    run.time = list(
      start = start.time,
      end = end.time,
      elapsed = difftime(end.time, start.time)
    ),
    reps = reps
  )
  saveRDS(res, res$filename)
  
  cat('\n--------', format(Sys.time()), 'End MAMBO --------\n')
  cat('  Number of replicates:', nrep, '\n')
  cat('  MCMC parameters:\n')
  cat('    Chains:', chains, '\n')
  cat('    Adapt:', adapt, '\n')
  cat('    Burnin:', burnin, '\n')
  cat('    Total Samples:', total.samples, '\n')
  cat('    Thinning:', thin, '\n')
  cat(
    '  Total elapsed time: ', 
    format(round(swfscMisc::autoUnits(res$run.time$elapsed), 1)),
    '\n',
    sep = ''
  )
  cat('  Results saved to:', res$filename, '\n')
  
  closeAllConnections()
  capture.output(gc(verbose = FALSE))
  invisible(res)
}
