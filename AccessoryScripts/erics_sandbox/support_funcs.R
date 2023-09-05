betaParams <- function(x) {
  # by-sample (columns) coverage
  coverage <- colSums(x)
  
  # fit beta shape parameters (transpose matrix for recycling of coverage vector)
  tx <- t(x)
  result <- array(c(tx + 1, coverage - tx + 1), dim = c(dim(tx), 2))
  dimnames(result) <- list(
    colnames(x),
    rownames(x),
    c('shape1', 'shape2')
  )
  result
}

ranRelPct <- function(beta.params) {
  # draw random sample from beta distribution with shape parameters 'p'
  pct <- apply(beta.params, c(1, 2), function(p) rbeta(1, p[1], p[2]))
  
  # normalize random percents to unity and return matrix
  # (transposed back to original dimensions)
  pct <- t(pct / rowSums(pct))
  rownames(pct) <- dimnames(beta.params)[[2]]
  colnames(pct) <- dimnames(beta.params)[[1]]
  pct
}

ranPCA <- function(beta.params) {
  prob <- t(ranRelPct(beta.params))
  log(prob / (1 - prob)) |> 
    prcomp() |> 
    summary()
}

numImpPCs <- function(pca) {
  # number of important PCs = as many as account for expected variance 
  imp.gt.exp <- pca$importance["Proportion of Variance", ] >= (1 / ncol(pca$importance))
  max(1, which.min(imp.gt.exp) - 1) 
}

jagsPClm <- function(pc.resp, pc.preds, chains, adapt, burnin, total.samples, thin) {
  library(runjags)
  
  data.list <- list(
    num.ind = nrow(pc.resp),
    num.resp = ncol(pc.resp),
    num.preds = ncol(pc.preds),
    pc.resp = pc.resp,
    pc.preds = pc.preds
  )
  
  jags.model <- "model {
    for(p in 1:num.resp) {
      # Prior for intercept
      intercept[p] ~ dnorm(0, 5e-4)

      for(a in 1:num.preds) {
        # Draw of switch for coefficient occurrence
        w[a, p] ~ dbern(0.5)
        
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
      for(p in 1:num.resp) {
        mu[i, p] <- intercept[p] + inprod(b.prime[, p], pc.preds[i, ])
        pc.resp[i, p] ~ dnorm(mu[i, p], tau[p])
      }
    }
  }"
  
  run.jags(
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

extractPC <- function(pca.list, pc, locus) {
  loadings <- sapply(pca.list, function(x) x$pca[[locus]]$rotation[, pc]) 
  for(i in 1:nrow(loadings)) {
    switch.sign <- sign(loadings[i, ]) != sign(loadings[i, 1])
    loadings[i, switch.sign] <- loadings[i, switch.sign] * -1
  }
  
  scores <- sapply(pca.list, function(x) x$pca[[locus]]$x[, pc])
  for(i in 1:nrow(scores)) {
    switch.sign <- sign(scores[i, ]) != sign(scores[i, 1])
    scores[i, switch.sign] <- scores[i, switch.sign] * -1
  }
  
  list(loadings = loadings, scores = scores)
}

#' Returns a vector of logicals identifying values that are 
#'   outliers. iqr = inter-quartile range, z = z-score
isOutlier <- function(x, type = c('iqr', 'z'), thresh = 3) {
  switch(
    match.arg(type),
    iqr = {
      quarts <- quantile(x, probs = c(0.25, 0.75))
      iqr <- diff(quarts)
      thresh <- c(quarts[1] - 1.5 * iqr, quarts[2] + 1.5 * iqr)
      x <= thresh[1] | x >= thresh[2]
    },
    z = {
      z.score <- (x - mean(x)) / sd(x)
      abs(z.score) >= thresh
    },
    NULL
  )
}

outlierLoadings <- function(pca) {
  apply(pca$rotation[, 1:pca$num.pcs], 2, function(x) {
    outliers <- x[isOutlier(x)]
    list(
      pos = sort(outliers[outliers > 0], decreasing = TRUE),
      neg = sort(outliers[outliers < 0]))
  })
}


contrastSummary <- function(results, d, pc, min.iter = length(results)) {
  res <- lapply(results, function(x) {
    outlierLoadings(x$pca[[d]])[[pc]]
  })
  
  pos <- lapply(res, function(x) names(x$pos)) |> 
    unlist() |> 
    table() |> 
    sort(d = T)
  
  neg <- lapply(res, function(x) names(x$neg)) |> 
    unlist() |> 
    table() |> 
    sort(d = T)
  
  list(
    pos = names(pos[pos >= min.iter]),
    neg = names(neg[neg >= min.iter])
  )
}


switchSummary <- function(results, min.p = 0.75) {
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
  
  smry <- w.post |> 
    as.data.frame.table(responseName = 'w') |> 
    mutate(rep = as.numeric(rep)) |> 
    group_by(pc.16s, pc.18s) |> 
    summarize(median = median(w), .groups = 'drop') |> 
    filter(median > min.p)
  
  p <- w.post |> 
    as.data.frame.table(responseName = 'w') |> 
    mutate(rep = as.numeric(rep)) |> 
    ggplot() +
    geom_histogram(aes(w)) +
    facet_grid(pc.16s ~ pc.18s)
  print(p)
  
  smry
}