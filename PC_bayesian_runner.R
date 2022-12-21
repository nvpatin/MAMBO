rm(list = ls())
source("ranRelPct.R")

# Load data ---------------------------------------------------------------
fl16s <- read.delim("Data/Flyer2018_16S_table_counts.tsv", row.names = 1)
fl18s <- read.delim("Data/Flyer2018_18S_table_counts.tsv", row.names = 1)

# Extract 16s PCs ---------------------------------------------------------
fl16s.prob <- t(ranRelPct(fl16s))
fl16s.lo <- log(fl16s.prob / (1 - fl16s.prob))
# using 'summary()' function to get variances of each component
fl16s.pca <- summary(prcomp(fl16s.lo))
fl16s.imp <- fl16s.pca$importance

# Extract 18s PCs ---------------------------------------------------------
# !!!! Make sure to sort the the 18s data so the samples are in the same order
# as the 16s data
fl18s.prob <- t(ranRelPct(fl18s))[rownames(fl16s.prob), ]
fl18s.lo <- log(fl18s.prob / (1 - fl18s.prob))
fl18s.pca <- summary(prcomp(fl18s.lo))

# These objects need to be in the workspace for the JAGS bayesian model to run
# number of samples
num.ind <- nrow(fl18s.pca$x)
# number of PCs to predict in 18S data
num.18s.pc <- 2
# number of predictor PCs to use from 16S data
# = as many as account for expected variance of predictor
num.preds <- which.min(fl16s.imp["Proportion of Variance", ] >= (1 / ncol(fl16s.imp)))
num.preds <- max(1, num.preds - 1)
# matrix of 18s PCs (rows are samples, columns are PCs)
pc.18s <- fl18s.pca$x[, 1:num.18s.pc]
# matrix of 16s PCs (rows are samples, columns are PCs)
pc.16s <- fl16s.pca$x[, 1:num.preds]

# MCMC parameters - need to be in workspace for model to run
chains <- 10
adapt <- 500
burnin <- 10000
total.samples <- 10000
thin <- 10

source("JAGS_PC_lm.R")

print(post$timetaken)
print(post.smry)
