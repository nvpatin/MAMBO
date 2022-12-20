#! /usr/local/bin/Rscript

rm(list = ls())
lasker2018 <- read.delim("Data/Lasker2018_table_counts.tsv", row.names = 1)

x <- lasker2018[1:40, 1]






#' Draws a random sample of the probability of occurrence of each OTU in each sample
sampleProbOccur <- function(x, num.cores = NULL) {
  #' fit dirichlet multinomial shape parameters
  #' (should be done once externally later)
  dm.params <- apply(x, 2, function(counts) {

  })

  # draw one random sample from dirichlet
  sample.dm <- function(alpha) rdirichlet(1, alpha)

  # draw matrix of random percent occurrence
  # use parallel processing if available
  if(is.null(num.cores)) num.cores <- parallel::detectCores() - 1
  if(num.cores > 1) {
    # spread draws across cores
    cl <- if(.Platform$OS.type == "windows") {
      parallel::makePSOCKcluster(num.cores)
    } else {
      parallel::makeForkCluster(num.cores)
    }
    prob.mat <- parallel::parApply(cl, dm.params, 2, sample.dm)
    parallel::stopCluster(cl)
    prob.mat
  } else {
    # single core mode
    apply(dm.params, 2, sample.dm)
  }
}

