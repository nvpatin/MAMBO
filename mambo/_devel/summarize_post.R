rm(list = ls())
library(mambo)

results <- readRDS('mambo.20231031_112830.rds')

w.smry <- switchSummary(results, 0.6)
w.smry$smry
