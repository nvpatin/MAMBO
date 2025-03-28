rm(list = ls())
library(mambo)
library(tidyverse)

fl16s <- read_tsv('../../Data/ASV tables/Flyer2018_16S_table_counts.tsv') |> 
  column_to_rownames('ASV ID') |> 
  filterCountData()

fl18s <- read_tsv('../../Data/ASV tables/Flyer2018_18S_table_counts.tsv') |> 
  column_to_rownames('ASV ID')|> 
  filterCountData()

result <- mambo(
  resp.label = '18s', 
  resp.counts = fl18s, 
  pred.label = '16s', 
  pred.counts = fl16s, 
  nrep = 100, 
  chains = 10
)

plotPCs(result, '16s')

plotPCs(result, '18s')
