rm(list = ls())
source('support_funcs.R')

fl16s <- betaParams('../../Data/Flyer2018_16S_table_counts.tsv')
fl18s <- betaParams('../../Data/Flyer2018_18S_table_counts.tsv')

res <- repBayesianPCAlm(
  run.label = 'function.test', 
  resp.label = '18s',
  resp.beta = fl18s,
  pred.label = '16s',
  pred.beta = fl16s, 
  nrep = 20, 
  mcmc = list(
    chains = 5, 
    adapt = 500,
    burnin = 1e4,
    total.samples = 5e3,
    thin = 100
  )
)tt