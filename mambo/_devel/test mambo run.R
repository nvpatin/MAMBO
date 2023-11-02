rm(list = ls())
library(mambo)

mambo(
  resp.label = '18s',
  resp.beta = betaParams(fl18s),
  pred.label = '16s',
  pred.beta = betaParams(fl16s),
  nrep = 10
)
