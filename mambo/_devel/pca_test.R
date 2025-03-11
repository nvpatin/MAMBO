rm(list = ls())
gc()
library(mambo)
library(tidyverse)

resp.counts <- read.csv('GOMECC Data/GOMECC_18S_ASV_table.csv') |> 
  column_to_rownames('hash') |> 
  filterCountData()

pred.counts <- read.csv('GOMECC Data/GOMECC_16S_ASV_table-V2.csv') |> 
  column_to_rownames('hash') |> 
  filterCountData()

x <- mambo(
  '18s', filterCountData(counts.18s), 
  '16s', filterCountData(counts.16s),
  nrep = 10, 
  bayesian = T
)

plotPCs(x, '16s')
