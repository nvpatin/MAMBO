#! /usr/local/bin/Rscript

rm(list = ls())
library(here)
source(here::here("ranRelPct.R"))

# For accepting input arguments
args <- commandArgs(trailingOnly=TRUE)

# Run with Lasker 2018 data --------------------------------------------------------
name <- args[3]
#name <- "Flyer2018_18S_table_counts.tsv"
path <- here::here("Data", name)
lasker2018 <- read.delim(path, row.names = 1)

ran.lasker <- ranRelPct(lasker2018)
colnames(ran.lasker) <- NULL
rownames(ran.lasker) <- NULL

options(max.print=2510000)
print(ran.lasker)
