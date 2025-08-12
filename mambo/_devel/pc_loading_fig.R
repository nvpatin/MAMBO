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

results <- readRDS('mambo.20250812_172628.rds')


plotLoadingBySource <- function(results, count.data, locus, pc, df) {
  ol <- outlierLoadings(results, locus, pc, thresh = 3.5)
  t(count.data)[, ol$asv] |> 
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
    labs(x = 'Sample', y = 'ASV', title = paste0(locus, ': PC', pc)) +
    facet_grid(sign ~ source, scales = 'free') +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 90, hjust = 0),
      legend.position = 'top',
      legend.key.width = unit(0.1, 'npc')
    ) 
}


pdf('pc loading summary.pdf', height = 12, width = 16)
plotPCs(results, '16s', sample.df = df, ellipse.fill = 'diversity', facet.by = 'source')
plotLoadingBySource(results, counts.16s, '16s', 1, df)
plotLoadingBySource(results, counts.16s, '16s', 2, df)
plotPCs(results, '18s', sample.df = df, ellipse.fill = 'diversity', facet.by = 'source')
plotLoadingBySource(results, counts.18s, '18s', 1, df)
plotLoadingBySource(results, counts.18s, '18s', 2, df)
dev.off()