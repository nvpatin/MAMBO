rm(list = ls())
library(mambo)

df <- data.frame(sample = colnames(counts.16s)) |> 
  mutate(
    cruise = case_when(
      stringr::str_starts(sample, 'Lasker') ~ 'Lasker',
      stringr::str_starts(sample, 'CN18S') ~ 'CN18S',
      stringr::str_starts(sample, 'CN18F') ~ 'CN18F'
    )
  ) 

by.cruise <- lapply(split(df, df$cruise), function(cruise.df) {
  mambo(
    resp.label = '18s', 
    resp.counts = counts.18s[, cruise.df$sample], 
    pred.label = '16s', 
    pred.counts = counts.16s[, cruise.df$sample], 
    nrep = 100, 
    bayesian = FALSE,
    run.label = unique(cruise.df$cruise)
  )
})

pdf('by_cruise_pca.pdf', height = 11, width = 11)
imap(by.cruise, function(x, idx) {
  div.df <- x$reps |> 
    lapply(function(x) x$sample.diversity) |> 
    dplyr::bind_rows() |> 
    pivot_longer(-sample, names_to = 'locus', values_to = 'diversity') |> 
    group_by(sample, locus) |> 
    summarize(diversity = median(diversity), .groups = 'drop') |> 
    pivot_wider(names_from = 'locus', values_from = 'diversity')
  for(locus in c('16s', '18s')) {
    p <- plotPCs(x, locus, sample.df = div.df, ellipse.fill = locus, plot = FALSE)
    print(p + ggtitle(paste0(idx, ': ', locus)))
  }
})
dev.off()
