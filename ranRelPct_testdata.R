#! /usr/local/bin/Rscript
library(here)

rm(list = ls())
source("/Users/nastassia.patin/GitHub/NOAA-NCAR-Hackathon/ranRelPct.R")

# Run with Lasker 2018 data --------------------------------------------------------
path <- here::here("Data", "Lasker2018_table_counts.tsv")
lasker2018 <- read.delim(path, row.names = 1)

ran.lasker <- ranRelPct(lasker2018)

options(max.print=2510000)
print(ran.lasker)
