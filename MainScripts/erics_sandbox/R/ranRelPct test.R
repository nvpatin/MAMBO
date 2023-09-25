rm(list = ls())
library(tidyverse)
library(runjags)

counts <- read.delim("../Data/Flyer2018_12S_table_counts.tsv") |>
  column_to_rownames("ASV.ID") |>
  as.matrix() |>
  t()

asvs.to.keep <- apply(counts, 2, function(x) any(x > 0))
mat <- counts[, asvs.to.keep]


mat <- matrix(sample(0:5, 42, T), nrow = 6)

chains <- 3
adapt <- 500
burnin <- 1e5
total.samples <- 1e3
thin <- 10

post <- run.jags(
  data = list(
    counts = mat,
    total = rowSums(mat),
    num.samples = nrow(mat),
    num.asvs = ncol(mat)
  ),
  model = "model{
    for(i in 1:num.samples) {
      for(j in 1:num.asvs) {
        alpha[i, j] ~ dunif(0, 40)
      }
      p[i, 1:num.asvs] ~ ddirch(alpha[i, 1:num.asvs])
      counts[i, ] ~ dmulti(p[i, 1:num.asvs], total[i])
    }
  }",
  monitor = c("deviance", "p", "alpha"),
  inits = function() list(
    .RNG.name = "lecuyer::RngStream",
    .RNG.seed = sample(1:9999, 1)
  ),
  modules = c("glm", "lecuyer"),
  summarise = FALSE,
  jags.refresh = 10,
  method = "parallel",
  n.chains = chains,
  adapt = adapt,
  burnin = burnin,
  sample = ceiling(total.samples / chains),
  thin = thin
)
end.time <- Sys.time()
elapsed <- swfscMisc::autoUnits(post$timetaken)

p <- swfscMisc::runjags2list(post)
