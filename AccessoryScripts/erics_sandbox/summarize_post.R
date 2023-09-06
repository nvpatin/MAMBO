rm(list = ls())
library(tidyverse)
source('support_funcs.R')

results <- readRDS('function.test_20230906_101154.rds')

w.smry <- switchSummary(results)
w.smry