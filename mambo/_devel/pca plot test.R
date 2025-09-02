rm(list = ls())
library(tidyverse)
library(mambo)

results <- readRDS('Lasker.20250819_124325.rds')

pca <- extractPCA(results)
x <- switchSummary(results)

plotPCs(results, 'Lasker_16S', 'Lasker_18S', 1, 1)
plotPCs(results, 'Lasker_16S', 'Lasker_18S', 2, 2)
plotPCs(results, 'Lasker_16S', 'Lasker_18S', 3, 3)
plotPCs(results, 'Lasker_16S', 'Lasker_18S', 4, 4)
plotPCs(results, 'Lasker_16S', 'Lasker_18S', 4, 6)
