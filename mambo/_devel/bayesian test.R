rm(list = ls())
library(mambo)
library(tidyverse)

num.cores <- parallel::detectCores() - 1

pred.counts <- read.csv('GOMECC Data/GOMECC_16S_ASV_table-V2.csv') |> 
  column_to_rownames('hash') |> 
  filterCountData()
pred.label <- '16s'

resp.counts <- read.csv('GOMECC Data/GOMECC_18S_ASV_table.csv') |> 
  column_to_rownames('hash')|> 
  filterCountData()
resp.label <- '18s'

result <- mambo(
  resp.label = '18s',
  resp.counts = resp.counts,
  pred.label = '16s',
  pred.counts = pred.counts,
  nrep=100,
  chains=3,
  adapt=100,
  burnin=1000,
  total.samples=1000,
  num.resp.pcs=5,
  num.pred.pcs=10,
  thin=1,
  run.label="mambo",
  output.log=TRUE
)

pcaSummary(result)
