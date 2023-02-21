rm(list = ls())
library(runjags)
source("ranRelPct.R")

# Load data ---------------------------------------------------------------
fl16s <- read.delim("../../../Data/Flyer2018_16S_table_counts.tsv", row.names = 1)
fl18s <- read.delim("../../../Data/Flyer2018_18S_table_counts.tsv", row.names = 1)

# Extract 16s PCs ---------------------------------------------------------
fl16s.prob <- t(ranRelPct(fl16s))
fl16s.lo <- log(fl16s.prob / (1 - fl16s.prob))
# using 'summary()' function to get variances of each component
fl16s.pca <- summary(prcomp(fl16s.lo))
fl16s.imp <- fl16s.pca$importance

# Extract 18s PCs ---------------------------------------------------------
# !!!! Make sure to sort the the 18s data so the samples are in the same order
# as the 16s data
fl18s.prob <- t(ranRelPct(fl18s))[rownames(fl16s.prob), ]
fl18s.lo <- log(fl18s.prob / (1 - fl18s.prob))
fl18s.pca <- summary(prcomp(fl18s.lo))

# number of predictor PCs to use from 16S data
# = as many as account for expected variance of predictor
num.preds <- which.min(fl16s.imp["Proportion of Variance", ] >= (1 / ncol(fl16s.imp)))
num.preds <- max(1, num.preds - 1)

# Number of 18s PCs to predict
num.18s.pc <- 2

chains <- 10
adapt <- 500
burnin <- 1e5
total.samples <- 1e4
thin <- 100

post <- run.jags(
  data = list(
    num.ind = nrow(fl18s.pca$x),
    num.18s.pc = num.18s.pc,
    num.preds = num.preds,
    pc.18s = fl18s.pca$x[, 1:num.18s.pc],
    pc.16s = fl16s.pca$x[, 1:num.preds]
  ),
  model = "model {
    for(p in 1:num.18s.pc) {
      # Prior for intercept
      intercept[p] ~ dnorm(0, 5e-4)

      for(a in 1:num.preds) {
        # Prior for coefficient occurrence
        pr.b[a, p] ~ dunif(0, 1)
        
        # Draw of switch for coefficient occurrence
        w[a, p] ~ dbern(pr.b[a, p])
        
        # Prior for ASV coefficients
        b[a, p] ~ dnorm(0, 1e-6)
        
        # Modified (switched) ASV coefficient
        b.prime[a, p] <- b[a, p] * w[a, p]
      }

      # Prior for variance
      v[p] ~ dunif(0, 1000)
      tau[p] <- 1 / v[p]
    }

    for(i in 1:num.ind) {
      for(p in 1:num.18s.pc) {
        mu[i, p] <- intercept[p] + inprod(b.prime[, p], pc.16s[i, ])
        pc.18s[i, p] ~ dnorm(mu[i, p], tau[p])
      }
    }
  }",
  monitor = c("deviance", "intercept", "pr.b", "b.prime",  "v"),
  inits = function() list(
    .RNG.name = "lecuyer::RngStream",
    .RNG.seed = sample(1:9999, 1)
  ),
  modules = c("glm", "lecuyer"),
  summarise = FALSE,
  method = "parallel",
  n.chains = chains,
  adapt = adapt,
  burnin = burnin,
  sample = ceiling(total.samples / chains),
  thin = thin
)
end.time <- Sys.time()
post$timetaken <- swfscMisc::autoUnits(post$timetaken)

print(post$timetaken)

plot(post, vars = "deviance")
