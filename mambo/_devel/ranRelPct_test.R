nrep <- 10
resp.counts <- asv.18s
pred.coutns <- asv.16s

# compute beta parameter arrays
resp.beta <- betaParams(resp.counts)
pred.beta <- betaParams(pred.counts)

# make sure rows are in same order for both sets of data
resp.beta <- resp.beta[dimnames(pred.beta)[[1]], , ]

reps <- lapply(1:nrep, function(i) {
  gc()
  start.time <- Sys.time()
  cat('\n-------- Replicate ')
  cat(i, '/', nrep, sep = '')
  cat(' --------\n')
  
  # Extract PCs -------------------------------------------------------------
  cat(' ', format(Sys.time()), 'ranRelPct...\n')
  
  resp.rel.pct <- ranRelPct(resp.beta)
  pred.rel.pct <- ranRelPct(pred.beta)
  
  end.time <- Sys.time()
  elapsed <- difftime(end.time, start.time)
  cat(
    '  End replicate:', 
    format(round(swfscMisc::autoUnits(elapsed))),
    '\n'
  )
  list(resp = resp.rel.pct, pred = pred.rel.pct)
})