#! /usr/local/bin/Rscript

rm(list = ls())
library(here)
source(here::here("ranRelPct.R"))

# For accepting input arguments
args <- commandArgs(trailingOnly=TRUE)

# Run with Lasker 2018 data --------------------------------------------------------
name <- args[3]
path <- here::here("Data", name)
lasker2018 <- read.delim(path, row.names = 1)

ran.lasker <- ranRelPct(lasker2018)

options(max.print=2510000)
print(ran.lasker)
