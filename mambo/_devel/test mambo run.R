rm(list = ls())
library(mambo)

lasker.samples <- sample.meta$sample[sample.meta$cruise == 'Lasker']

lasker.16s <- filterCountData(counts.16s[, lasker.samples])
lasker.18s <- filterCountData(counts.18s[, lasker.samples])

mambo(
  resp.label = 'Lasker_18S', 
  resp.counts = lasker.18s, 
  pred.label = 'Lasker_16S', 
  pred.counts = lasker.16s, 
  nrep = 10, 
  chains = 3,
  adapt = 100,
  burnin = 1000,
  run.label = 'Lasker'
)
