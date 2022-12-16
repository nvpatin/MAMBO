# Function to draw 'n' random relative percent matrices -------------------

ranRelPct <- function(n, num.reads) {
  coverage <- matrix(
    rep(colSums(num.reads), each = nrow(num.reads)),
    ncol = ncol(num.reads)
  )

  # fit beta shape parameters
  beta.params <- array(
    c(num.reads + 1, coverage - num.reads + 1),
    dim = c(dim(num.reads), 2)
  )

  # draw n random matrices from beta distribution
  replicate(n, {
    pct <- apply(beta.params, c(1, 2), function(x) rbeta(1, x[1], x[2]))
    t(t(pct) / colSums(pct))
  }, simplify = FALSE)
}
