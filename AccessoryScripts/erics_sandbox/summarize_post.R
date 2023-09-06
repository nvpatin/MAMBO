rm(list = ls())
library(tidyverse)
source('support_funcs.R')

results <- readRDS('function.test_20230905_190026.rds')

w.smry <- switchSummary(results)
w.smry