rm(list = ls())
library(tidyverse)

results <- readRDS('bayesian.reps.rds')

min.pcs <- results |> 
  sapply(function(x) {
    sapply(x$pca, function(pca.x) pca.x$num.pcs)
  }) |> 
  apply(1, min)

w.post <- do.call(
  abind::abind,
  c(lapply(results, function(x) {
    apply(x$post.list$w[1:min.pcs['16s'], 1:min.pcs['18s'], ], c(1, 2), mean)
  }), list(along = 3))
)
dimnames(w.post)[[3]] <- 1:dim(w.post)[3]
names(dimnames(w.post)) <- c('pc.16s', 'pc.18s', 'rep')
  
w.post |> 
  as.data.frame.table(responseName = 'w') |> 
  mutate(rep = as.numeric(rep)) |> 
  ggplot() +
  geom_histogram(aes(w)) +
  facet_grid(pc.16s ~ pc.18s)

w.post |> 
  as.data.frame.table(responseName = 'w') |> 
  mutate(rep = as.numeric(rep)) |> 
  group_by(pc.16s, pc.18s) |> 
  summarize(median = median(w), .groups = 'drop') |> 
  filter(median > 0.7)

