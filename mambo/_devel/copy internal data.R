rm(list = ls())
library(tidyverse)

counts.16s <- read.csv(
  '../../Data/ASV tables/CCE/Merged2018_16S_otu_filtered.csv',
  row.names = 1
)
save(counts.16s, file = '../data/counts.16s.rda')

counts.18s <- read.csv(
  '../../Data/ASV tables/CCE/Merged2018_18S_otu_filtered.csv',
  row.names = 1
)
save(counts.18s, file = '../data/counts.18s.rda')

counts.12s <- read.csv(
  '../../Data/ASV tables/Blackman/p571_run200218_12S_ASVs.csv',
  row.names = 1
)
save(counts.12s, file = '../data/counts.12s.rda')

counts.COI <- read.csv(
  '../../Data/ASV tables/Blackman/p571_run200121_COI_ASVs.csv',
  row.names = 1
)
save(counts.COI, file = '../data/counts.COI.rda')

taxa.16s <- read.csv(
  '../../Data/Taxonomy tables/CCE/Merged2018_16S_taxa_filtered.csv',
  row.names = 1
)
save(taxa.16s, file = '../data/taxa.16s.rda')

taxa.18s <- read.csv(
  '../../Data/Taxonomy tables/CCE/Merged2018_18S_taxa_filtered.csv',
  row.names = 1
)
save(taxa.18s, file = '../data/taxa.18s.rda')

taxa.12s <- read.csv(
  '../../Data/Taxonomy tables/Blackman/p571_run200218_12S_taxonomy.csv',
  row.names = 1
)
save(taxa.12s, file = '../data/taxa.12s.rda')

taxa.COI <- read.csv(
  '../../Data/Taxonomy tables/Blackman/p571_run200121_COI_taxonomy.csv',
  row.names = 1
)
save(taxa.COI, file = '../data/taxa.COI.rda')

sample.meta <- read.csv('../../Data/Metadata/Merged2018_meta_filtered.csv') |> 
  rename(sample = 'sample_name')
save(sample.meta, file = '../data/sample.meta.rda')
