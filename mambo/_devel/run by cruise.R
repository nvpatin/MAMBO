rm(list = ls())
library(mambo)

df <- data.frame(sample = colnames(counts.16s)) |> 
  mutate(
    source = case_when(
      stringr::str_starts(sample, 'Lasker') ~ 'Lasker',
      stringr::str_starts(sample, 'CN18S') ~ 'CN18S',
      stringr::str_starts(sample, 'CN18F') ~ 'CN18F'
    )
  ) 

by.source <- lapply(split(df, df$source), function(source.df) {
  mambo(
    resp.label = '18s', 
    resp.counts = counts.18s[, source.df$sample], 
    pred.label = '16s', 
    pred.counts = counts.16s[, source.df$sample], 
    nrep = 50, 
    bayesian = FALSE,
    run.label = unique(source.df$source)
  )
})

pdf('by_cruise_pca.pdf', height = 11, width = 11)
imap(by.source, function(x, idx) {
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
