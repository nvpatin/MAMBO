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
  resp.label = resp.label,
  resp.counts = resp.counts,
  pred.label = pred.label,
  pred.counts = pred.counts,
  num.resp.pcs = 5,
  num.pred.pcs = 5,
  nrep = 10,
  chains = 10
)

x <- pcaSummary(result)