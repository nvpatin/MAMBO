rm(list = ls())
library(mambo)
library(tidyverse)

results <- readRDS('Lasker.20250819_124325.rds')

plotPCs2(result, 'Lasker_18S', sample.df = sample.meta, ellipse.fill = 'depth_m')
