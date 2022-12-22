#! /usr/local/bin/Rscript

rm(list = ls())
library(here)
source(here::here("ranRelPct.R"))

# For accepting input arguments
args <- commandArgs(trailingOnly=TRUE)

# Run with Flyer 2018 data --------------------------------------------------------
name <- args[3]
path <- here::here("Data", name)
lasker2018 <- read.delim(path, row.names = 1)

# Sort columns alphabetically so 16S and 18S data frames will be in the same order
new_order = sort(colnames(lasker2018))
lasker2018 <- lasker2018[, new_order]

ran.lasker <- ranRelPct(lasker2018)

# Remove column and row names
colnames(ran.lasker) <- NULL
rownames(ran.lasker) <- NULL

# Transpose matrix for Python import
ran.lasker.forpy <- as.vector(t(ran.lasker))

options(max.print=500000)
cat(ran.lasker.forpy)
