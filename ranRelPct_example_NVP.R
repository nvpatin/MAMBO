#! /usr/local/bin/Rscript
library(here)

rm(list = ls())
source("/Users/nastassia.patin/GitHub/NOAA-NCAR-Hackathon/ranRelPct.R")

# Import ASV counts sheet ----------------------------------------

path <- here::here("GitHub", "NOAA-NCAR-Hackathon", 
                        "Data", "Lasker2018_table_counts.tsv")
num.reads.2 <- as.matrix(read.csv(path, row.names = 1, sep='\t'))

# Test it --------------------------------------------------------

# This is the mode (expected value) of the distribution
rel.pct <- t(t(num.reads.2) / colSums(num.reads.2))

# One random draw
ran.1 <- ranRelPct(1, num.reads.2)

# Ten random draws
ran.10 <- ranRelPct(10, num.reads)

print(ran.1, sep=',')
