rm(list = ls())
library(mambo)

results <- readRDS('mambo.20231031_112830.rds')

ol1.df <- outlierLoadings(results, '18s', 1)
head(ol1.df, 10)
tail(ol1.df, 10)

w.smry <- switchSummary(results, 0.6)
w.smry$smry

plotPCs(results, '18s')