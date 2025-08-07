rm(list = ls())
library(mambo)
library(tidyverse)

result <- mambo(
  resp.label = '18s', 
  resp.counts = counts.18s, 
  pred.label = '16s', 
  pred.counts = counts.16s, 
  nrep = 10, 
  bayesian = FALSE
)
