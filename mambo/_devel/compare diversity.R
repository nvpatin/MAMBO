rm(list = ls())
library(tidyverse)
library(mambo)

results <- readRDS('mambo.20250812_172628.rds')

df <- data.frame(sample = colnames(counts.16s)) |> 
  mutate(
    source = case_when(
      stringr::str_starts(sample, 'Lasker') ~ 'Lasker',
      stringr::str_starts(sample, 'CN18S') ~ 'CN18S',
      stringr::str_starts(sample, 'CN18F') ~ 'CN18F'
    )
  ) 

pdf('diversity summary.pdf', height = 11, width = 11)
results$reps |> 
  lapply(function(x) x$total.diversity) |> 
  dplyr::bind_rows() |> 
  pivot_longer(everything(), names_to = 'locus', values_to = 'diversity') |> 
  ggplot() +
  geom_histogram(aes(diversity), alpha = 0.7, bins = 20) +
  scale_x_log10() +
  facet_grid(~ locus, scales = 'free') 

results$reps |> 
  lapply(function(x) x$sample.diversity) |> 
  dplyr::bind_rows() |> 
  mutate(
    source = case_when(
      stringr::str_starts(sample, 'Lasker') ~ 'Lasker',
      stringr::str_starts(sample, 'CN18S') ~ 'CN18S',
      stringr::str_starts(sample, 'CN18F') ~ 'CN18F'
    )
  ) |> 
  pivot_longer(-c(sample, source), names_to = 'locus', values_to = 'diversity') |> 
  ggplot() +
  geom_histogram(aes(diversity, fill = source), alpha = 0.7, bins = 20) +
  scale_x_log10() +
  facet_grid(source ~ locus, scales = 'free') +
  theme(legend.position = 'top')

plotDiversity(results, sample.df = df, facet.by = 'source')
dev.off()