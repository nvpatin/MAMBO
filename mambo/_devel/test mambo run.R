rm(list = ls())
library(mambo)

result <- mambo(
  resp.label = '18s', 
  resp.counts = counts.18s, 
  pred.label = '16s', 
  pred.counts = counts.16s, 
  nrep = 5, 
  chains = 10
)
