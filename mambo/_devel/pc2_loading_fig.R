rm(list = ls())
library(tidyverse)
library(mambo)

df <- counts.16s |> 
  rownames_to_column('asv') |> 
  pivot_longer(-asv, names_to = 'sample', values_to = 'count') |> 
  mutate(
    source = case_when(
      stringr::str_starts(sample, 'Lasker') ~ 'Lasker',
      stringr::str_starts(sample, 'CN18S') ~ 'CN18S',
      stringr::str_starts(sample, 'CN18F') ~ 'CN18F'
    )
  ) |> 
  group_by(sample, source) |> 
  mutate(pct = count / sum(count)) |> 
  summarize(
    read.count = sum(count),
    diversity = sprex::diversity(pct),
    .groups = 'drop'
  )

results <- readRDS('mambo.20250809_230611.rds')

plotPCs(results, '16s', sample.df = df, facet.by = 'source')


ol <- outlierLoadings(results, '16s', 2)

pca <- extractPCA(results)

outlier.16s <- pca$scores[['16s']] |> 
  filter(pc == 1 & score > 60) |> 
  pull('sample') |> 
  unique()

cn18f <- df |> 
  filter(source == 'CN18F') |> 
  pull('sample')

cn18f.pc2 <- t(counts.16s)[cn18f, ol$asv]

cn18f.pc2 |> 
  as.data.frame() |> 
  rownames_to_column('sample') |> 
  pivot_longer(-sample, names_to = 'asv', values_to = 'counts') |>
  left_join(ol, by = 'asv') |> 
  left_join(df, by = 'sample') |> 
  mutate(
    sample = reorder(sample, diversity),
    asv = reorder(asv, abs(median.loading)),
    counts = ifelse(counts == 0, NA, counts),
    sign = ifelse(median.loading > 0, 'positive', 'negative'),
    sign = factor(sign, levels = c('positive', 'negative'))
  ) |> 
  ggplot(aes(sample, asv)) +
  geom_tile(aes(fill = counts)) +
  scale_fill_viridis_c() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0),
    legend.position = 'top',
    legend.key.width = unit(0.2, 'npc')
  ) +
  facet_wrap(~sign, ncol = 1, scales = 'free_y')
  

