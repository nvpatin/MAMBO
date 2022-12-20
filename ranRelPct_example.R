#! /usr/local/bin/Rscript
library(here)

rm(list = ls())
source("/Users/nastassia.patin/GitHub/NOAA-NCAR-Hackathon/ranRelPct.R")

# Example occurrence data -------------------------------------------------

#num.samples <- 5

#num.otus <- 10

#num.reads <- matrix(
#  sample(0:30, num.samples * num.otus, replace = TRUE),
#  ncol = num.samples
#)


# Test it -----------------------------------------------------------------

# This is the mode (expected value) of the distribution
#rel.pct <- t(t(num.reads) / colSums(num.reads))

# One random draw
#ran.pct <- ranRelPct(num.reads)


# Lasker 2018 data --------------------------------------------------------
path <- here::here("Data", "Lasker2018_table_counts.tsv")
lasker2018 <- read.delim(path, row.names = 1)

ran.lasker <- ranRelPct(lasker2018)

options(max.print=2510000)
print(ran.lasker)
