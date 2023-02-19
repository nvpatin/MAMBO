rm(list = ls())
library(runjags)

fl16s <- read.delim("../Data/Flyer2018_16S_table_counts.tsv", row.names = 1)

post <- run.jags(
  model = "model {
    for(i in 1:n.ind) {
      p[1:n.asv, i] ~ ddirch(alpha[1:n.asv])
      counts[, i] ~ dmulti(p[, i], total.reads[i])
    }
  }",
  monitor = c("deviance", "p"),
  data = list(
    n.ind = ncol(fl16s),
    n.asv = nrow(fl16s),
    counts = as.matrix(fl16s),
    total.reads = colSums(fl16s),
    alpha = rep(1, nrow(fl16s))
  ),
  inits = function() list(
    .RNG.name = "lecuyer::RngStream",
    .RNG.seed = sample(1:9999, 1)
  ),
  modules = c("glm", "lecuyer"),
  summarise = FALSE,
  jags.refresh = 10,
  method = "parallel",
  n.chains = 4,
  adapt = 1000,
  burnin = 10000,
  sample = 250,
  thin = 10
)

p <- myFuncs::runjags2list(post)
save(post, p, file = "count_multinomial_posterior.rdata")
