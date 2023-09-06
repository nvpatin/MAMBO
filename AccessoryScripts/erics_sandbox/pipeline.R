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


repBayesianPCAlm(
  run.label = 'function.test', 
  resp.label = '18s',
  resp.beta = fl18s,
  pred.label = '16s',
  pred.beta = fl16s, 
  nrep = 20, 
  mcmc = list(
    chains = 6, 
    adapt = 500,
    burnin = 1e4,
    total.samples = 5e3,
    thin = 50
  )
)