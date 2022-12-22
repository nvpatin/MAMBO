#! /usr/local/bin/Rscript

rm(list = ls())
library(here)
library(dplyr)
source(here::here("ranRelPct.R"))
source(here::here("jagsPClm.R"))

# For accepting input arguments
args <- commandArgs(trailingOnly=TRUE)

# Load data ---------------------------------------------------------------
fl16s.prob <- read.delim(here::here("Data", args[3]), row.names = 1)
fl18s.prob <- read.delim(here::here("Data", args[4]), row.names = 1)

# Extract 16s PCs ---------------------------------------------------------
fl16s.lo <- log(fl16s.prob / (1 - fl16s.prob))
# using 'summary()' function to get variances of each component
fl16s.pca <- summary(prcomp(fl16s.lo))
fl16s.imp <- fl16s.pca$importance

# Extract 18s PCs ---------------------------------------------------------
# !!!! Make sure to sort the the 18s data so the samples are in the same order
# as the 16s data
fl18s.lo <- log(fl18s.prob / (1 - fl18s.prob))
fl18s.pca <- summary(prcomp(fl18s.lo))

# number of PCs to predict in 18S data
num.18s.pc <- 2
# number of predictor PCs to use from 16S data
# = as many as account for expected variance of predictor
num.preds <- which.min(fl16s.imp["Proportion of Variance", ] >= (1 / ncol(fl16s.imp)))
num.preds <- max(1, num.preds - 1)

invisible(capture.output(post.smry <- jagsPClm(
  num.ind = nrow(fl18s.pca$x),
  num.18s.pc = num.18s.pc,
  num.preds = num.preds,
  pc.18s = fl18s.pca$x[, 1:num.18s.pc],
  pc.16s = fl16s.pca$x[, 1:num.preds],
  chains = 10,
  adapt = 500,
  burnin = 10000,
  total.samples = 10000,
  thin = 10
)))

df <- as.data.frame(post.smry)
df$rownames <- rownames(df)
df <- df %>% relocate(rownames, .before = Lower95)
df[46,] <- colnames(df)
df <- df[order(df$rownames=='rownames', decreasing=TRUE), ]

smry.forpy <- as.vector(t(df))
print(smry.forpy)
