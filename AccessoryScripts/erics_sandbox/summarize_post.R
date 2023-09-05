rm(list = ls())
library(tidyverse)
source('support_funcs.R')

results <- readRDS('bayesian.reps.rds')

w.smry <- switchSummary(results)
w.smry
