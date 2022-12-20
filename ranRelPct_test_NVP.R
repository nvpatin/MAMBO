#! /usr/local/bin/Rscript

# A test script to see if we can import a raw counts table, run \
# an operation (e.g. convert to relative abundances), and export 
# it back to Python

library(here)

rm(list = ls())
source("/Users/nastassia.patin/GitHub/NOAA-NCAR-Hackathon/ranRelPct.R")

# Example occurrence data -------------------------------------------------

num.samples <- 5
num.otus <- 10

num.reads <- matrix(
  sample(0:30, num.samples * num.otus, replace = TRUE),
  ncol = num.samples
)

#print(num.reads, sep=',')

# Import ASV counts sheet ----------------------------------------

path <- here::here("Data", "Lasker2018_table_counts.tsv")
num.reads.2 <- as.matrix(read.csv(path, row.names = 1, sep='\t'))

# This is the mode (expected value) of the distribution
rel.pct <- t(t(num.reads.2) / colSums(num.reads.2))

# One random draw
ran.1 <- ranRelPct(1, num.reads.2)

print(ran.1, sep=',')
