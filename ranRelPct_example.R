#! /usr/local/bin/Rscript
rm(list = ls())
source("ranRelPct.R")

# Example occurrence data -------------------------------------------------

num.samples <- 5

num.otus <- 10

num.reads <- matrix(
  sample(0:30, num.samples * num.otus, replace = TRUE),
  ncol = num.samples
)


# Test it -----------------------------------------------------------------

# This is the mode (expected value) of the distribution
rel.pct <- t(t(num.reads) / colSums(num.reads))

# One random draw
ran.pct <- ranRelPct(num.reads)

rel.pct
ran.pct


# Lasker 2018 data --------------------------------------------------------
lasker2018 <- read.delim("Data/Lasker2018_table_counts.tsv", row.names = 1)

ran.lasker <- ranRelPct(lasker2018)

# options(max.print=2510000)
head(ran.lasker, 20)
