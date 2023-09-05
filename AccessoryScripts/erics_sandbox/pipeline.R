rm(list = ls())
source('support_funcs.R')


# Load data ---------------------------------------------------------------
fl16s <- '../../Data/Flyer2018_16S_table_counts.tsv' |> 
  read.delim(row.names = 1) |> 
  betaParams()
fl18s <- '../../Data/Flyer2018_18S_table_counts.tsv' |> 
  read.delim(row.names = 1) |> 
  betaParams()
# make sure rows are in same order for both sets of data
fl18s <- fl18s[dimnames(fl16s)[[1]], , ]


# Multiple draws and Bayesian models --------------------------------------
result <- lapply(1:5, function(i) {
  cat('\n-------------\n')
  cat('Replicate', i, '\n')
  cat('-------------\n\n')
  
  # Extract PCs -------------------------------------------------------------
  pca <- list('16s' = ranPCA(fl16s), '18s' = ranPCA(fl18s))
  pca$'16s'$num.pcs = numImpPCs(pca$'16s')
  pca$'18s'$num.pcs = numImpPCs(pca$'18s')
  
  # Run Bayesian model ------------------------------------------------------
  post <- jagsPClm(
    pc.resp = pca$'18s'$x[, 1:pca$'18s'$num.pcs],
    pc.preds = pca$'16s'$x[, 1:pca$'16s'$num.pcs],
    chains = 6, 
    adapt = 500, 
    burnin = 3e4, 
    total.samples = 1e3, 
    thin = 10
  )
  
  cat('\n-------------\n')
  cat('Time elapsed:', swfscMisc::autoUnits(post$timetaken), '\n')
  cat('-------------\n\n')
  
  p <- swfscMisc::runjags2list(post)
  dimnames(p$intercept)[[1]] <-
    dimnames(p$b.prime)[[2]] <-
    dimnames(p$w)[[2]] <-
    dimnames(p$v)[[1]] <- paste0('18s.', 'PC', 1:pca$'18s'$num.pcs)
  dimnames(p$b.prime)[[1]] <-
    dimnames(p$w)[[1]] <- paste0('16s.', 'PC', 1:pca$'16s'$num.pcs)
  
  list(pca = pca, post.smry = summary(post), post.list = p)
})

saveRDS(result, 'bayesian.reps.rds')