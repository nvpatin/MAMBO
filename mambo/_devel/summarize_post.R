rm(list = ls())
library(mambo)

results <- readRDS('function.test_20230906_115436.rds')

w.smry <- switchSummary(results, 0.6)
w.smry
