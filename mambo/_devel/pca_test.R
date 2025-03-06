library(mambo)
load('../data/counts.18s.rda')
load('../data/counts.16s.rda')

nrep <- 4
resp.label <- '18s'
resp.counts <- counts.18s
pred.label <- '16s'
pred.counts <- counts.16s

# compute beta parameter arrays
resp.beta <- betaParams(resp.counts)
pred.beta <- betaParams(pred.counts)

# make sure rows are in same order for both sets of data
resp.beta <- resp.beta[dimnames(pred.beta)[[1]], , ]

reps <- lapply(1:nrep, function(i) {
  start.time <- Sys.time()
  cat('\n-------- Replicate ')
  cat(i, '/', nrep, sep = '')
  cat(' --------\n')
  
  # Extract PCs -------------------------------------------------------------
  cat(' ', format(Sys.time()), 'PCA...\n')
  pca <- stats::setNames(
    list(ranPCA(resp.beta), ranPCA(pred.beta)),
    c(resp.label, pred.label)
  )
  
  n <- numImpPCs(pca[[resp.label]])
  pca[[resp.label]]$rotation <- pca[[resp.label]]$rotation[, 1:n]
  pca[[resp.label]]$x<- pca[[resp.label]]$x[1:n, 1:n]
  pca[[resp.label]]$importance <- pca[[resp.label]]$importance[, 1:n]
  pca[[resp.label]]$num.pcs <- n
  pca[[resp.label]]$sdev <- pca[[resp.label]]$center <- pca[[resp.label]]$scale <- NULL
  
  n <- numImpPCs(pca[[pred.label]])
  pca[[pred.label]]$rotation <- pca[[pred.label]]$rotation[, 1:n]
  pca[[pred.label]]$x<- pca[[pred.label]]$x[1:n, 1:n]
  pca[[pred.label]]$importance <- pca[[pred.label]]$importance[, 1:n]
  pca[[pred.label]]$num.pcs <- n
  pca[[pred.label]]$sdev <- pca[[pred.label]]$center <- pca[[pred.label]]$scale <- NULL
  
  end.time <- Sys.time()
  elapsed <- difftime(end.time, start.time)
  cat(
    '  End replicate:', 
    format(round(swfscMisc::autoUnits(elapsed))),
    '\n'
  )
  pca
})
