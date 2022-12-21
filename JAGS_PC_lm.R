library(runjags)

post <- run.jags(
  model = "model {
    for(p in 1:num.18s.pc) {
      # Prior for intercept
      intercept[p] ~ dnorm(0, 5e-4)

      # Prior for ASV coefficients
      for(a in 1:num.preds) {
        b[a, p] ~ dnorm(0, 1e-6)
      }

      # Prior for variance
      v[p] ~ dunif(0, 1000)
      tau[p] <- 1 / v[p]
    }

    for(i in 1:num.ind) {
      for(p in 1:num.18s.pc) {
        mu[i, p] <- intercept[p] + inprod(b[, p], pc.16s[i, ])
        pc.18s[i, p] ~ dnorm(mu[i, p], tau[p])
      }
    }
  }",
  monitor = c("deviance", "intercept", "b", "v"),
  data = list(
    num.ind = num.ind,
    num.18s.pc = num.18s.pc,
    num.preds = num.preds,
    pc.18s = pc.18s,
    pc.16s = pc.16s
  ),
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

post.smry <- summary(post)
