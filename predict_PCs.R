rm(list = ls())
library(runjags)
source("ranRelPct.R")

fl16s <- read.delim("Data/Flyer2018_16S_table_counts.tsv", row.names = 1)
fl18s <- read.delim("Data/Flyer2018_18S_table_counts.tsv", row.names = 1)

fl16s.prob <- t(ranRelPct(fl16s))
fl16s.lo <- log(fl16s.prob / (1 - fl16s.prob))
fl18s.prob <- t(ranRelPct(fl18s))
fl18s.lo <- log(fl18s.prob / (1 - fl18s.prob))[rownames(fl16s.lo), ]

# rm(fl16s, fl18s, fl16s.prob, fl18s.prob, ranRelPct)
# gc()

fl16s.pca <- prcomp(fl16s.lo)
fl18s.pca <- prcomp(fl18s.lo)

model.data <- list(
  num.ind = nrow(fl18s.pca$x),
  num.pc = 2, #ncol(fl18s.pc),
  num.asv = ncol(fl16s.pca$x),
  fl18s.pc = fl18s.pca$x[, 1:2],
  fl16s.pc = fl16s.pca$x
)

post <- run.jags(
  model = "model {
    for(p in 1:num.pc) {
      # Prior for intercept
      intercept[p] ~ dnorm(0, 5e-4)

      # Prior for ASV coefficients
      for(a in 1:num.asv) {
        b[a, p] ~ dnorm(0, 1e-6)
      }

      # Prior for variance
      v[p] ~ dunif(0, 1000)
      tau[p] <- 1 / v[p]
    }

    for(i in 1:num.ind) {
      for(p in 1:num.pc) {
        mu[i, p] <- intercept[p] + inprod(b[, p], fl16s.pc[i, ])
        fl18s.pc[i, p] ~ dnorm(mu[i, p], tau[p])
      }
    }
  }",
  monitor = c("deviance", "intercept", "b", "v"),
  data = model.data,
  inits = function() list(
    .RNG.name = "lecuyer::RngStream",
    .RNG.seed = sample(1:9999, 1)
  ),
  modules = c("glm", "lecuyer"),
  summarise = FALSE,
  method = "parallel",
  n.chains = 8,
  adapt = 100,
  burnin = 100000,
  sample = 10000,
  thin = 1
)
end.time <- Sys.time()
post$timetaken <- swfscMisc::autoUnits(post$timetaken)

save.image("18s PC bayesian prediction.rdata")

graphics.off()
pdf("posterior plots.pdf")
plot(post, vars = c("deviance", "b"))
dev.off()
