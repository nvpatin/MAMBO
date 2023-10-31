#' @title Run JAGS linear models of principal components
#' @description Run a bayesian linear model with principal component response
#'   and predictors in JAGS.
#'
#' @param pc.resp matrix of response component scores.
#' @param pc.preds matrix of predictor component scores.
#' @param chains number of MCMC chains.
#' @param adapt number of adaptation iterations.
#' @param burnin number of burnin iterations.
#' @param total.samples total number of samples from the posterior to save.
#' @param thin number of iterations to skip between samples in each chain.
#'
#' @return a run.jags posterior mcmc object.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#'
jagsPClm <- function(pc.resp, pc.preds, chains, adapt, burnin, total.samples, thin) {
  data.list <- list(
    num.ind = nrow(pc.resp),
    num.resp = ncol(pc.resp),
    num.preds = ncol(pc.preds),
    pc.resp = pc.resp,
    pc.preds = pc.preds
  )

  jags.model <- "model {
    for(r in 1:num.resp) {
      # Prior for intercept
      intercept[r] ~ dnorm(0, 5e-4)

      for(p in 1:num.preds) {
        # Draw of switch for coefficient occurrence
        w[r, p] ~ dbern(0.5)

        # Prior for ASV coefficients
        b[r, p] ~ dnorm(0, 1e-6)

        # Modified (switched) ASV coefficient
        b.prime[r, p] <- b[r, p] * w[r, p]
      }

      # Prior for variance
      v[r] ~ dunif(0, 1000)
      tau[r] <- 1 / v[r]
    }

    for(i in 1:num.ind) {
      for(r in 1:num.resp) {
        mu[i, r] <- intercept[r] + inprod(b.prime[r, ], pc.preds[i, ])
        pc.resp[i, r] ~ dnorm(mu[i, r], tau[r])
      }
    }
  }"

  runjags::run.jags(
    model = jags.model,
    monitor = c("deviance", "intercept", "b.prime", "w", "v"),
    data = data.list,
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
}
