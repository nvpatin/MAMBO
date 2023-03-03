#! /usr/local/bin/Rscript

rm(list = ls())
library(here)
library(dplyr)
source(here::here("ranRelPct.R"))
source(here::here("jagsPClm.R"))

# For accepting input arguments
args <- commandArgs(trailingOnly=TRUE)

# Load data ---------------------------------------------------------------
pc.16s <- as.matrix(read.delim(here::here("Data", 
                                    "Flyer2018_16S_PCs.tsv"), 
                         row.names = 1))
pc.18s <- as.matrix(read.delim(here::here("Data", 
                                    "Flyer2018_18S_PCs.tsv"), 
                         row.names = 1))

# number of PCs to predict in 18S data
num.18s.pc <- strtoi(args[3])
# number of predictor PCs to use from 16S data
# = as many as account for expected variance of predictor
num.preds <- strtoi(args[4])
# number of observances (samples)
num.ind <- strtoi(args[5])

invisible(capture.output(post.smry <- jagsPClm(
  num.ind = num.ind,
  num.18s.pc = num.18s.pc,
  num.preds = num.preds,
  pc.18s = pc.18s,
  pc.16s = pc.16s,
  chains = 10,
  adapt = 500,
  burnin = 10000,
  total.samples = 10000,
  thin = 10
)))

df <- as.data.frame(post.smry)
df$rownames <- rownames(df)
df <- df %>% relocate(rownames, .before = Lower95)
df[18,] <- colnames(df) # This needs to change based on number of 16S predictors; 18 is for 6 predictors
df <- df[order(df$rownames=='rownames', decreasing=TRUE), ]

smry.forpy <- as.vector(t(df))
print(smry.forpy)
